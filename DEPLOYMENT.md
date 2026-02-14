# Deploying R Shiny app to AWS Lightsail/Container

2026-02-13

## Background

We've been having persistent issues with app stability on shinyapps.io hosting, so trying AWS Lightsail+Container, so straight deployment from local to AWS, and it takes care of build, etc.

Unless otherwise noted, all the following are local.

## Prepare the shiny app

### Create Dockerfile

Create a `Dockerfile` in the project root - see Dockerfile itself for details

## Build docker container

See `Dockerfile`, it contains the app and R packages, etc.

```sh
# Build the Docker image
docker build -t islets-shiny .

# Test locally, Visit http://localhost:3838
docker run -p 3838:3838 islets-shiny
```

## Configure AWS CLI
```sh
Enter your:

AWS Access Key ID
AWS Secret Access Key
Default region (e.g., us-east-1)
Default output format (e.g., json)
```

## Create Lightsail Container Service

```sh
# Create service with medium power (2GB RAM, 1 vCPU) - $40/month
aws lightsail create-container-service \
  --service-name islets-app \
  --power medium \
  --scale 1

# Check status (wait until state is "ACTIVE")
aws lightsail get-container-services --service-name islets-app \
  --query 'containerServices[0].state' \
  --output text
```

Note: Service provisioning takes 5-10 minutes. Wait until status shows ACTIVE before proceeding.

##  Push Docker Image to Lightsail

The standard Lightsail push command can have issues on WSL2, so we use Amazon ECR (Elastic Container Registry) instead.

```sh
# Get your AWS account ID
AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
echo ${AWS_ACCOUNT_ID}

# Create ECR repository
aws ecr create-repository --repository-name islets-shiny --region us-east-1

# Login to ECR
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com

# Tag your image for ECR
docker tag islets-shiny:latest ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/islets-shiny:latest

# Push to ECR (this may take 10-20 minutes for a 2.5GB image)
docker push ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/islets-shiny:latest

# Enable private registry access for Lightsail
aws lightsail update-container-service \
  --service-name islets-app \
  --private-registry-access '{"ecrImagePullerRole":{"isActive":true}}'

# Get the Lightsail ECR puller role ARN
LIGHTSAIL_ROLE=$(aws lightsail get-container-services \
  --service-name islets-app \
  --query 'containerServices[0].privateRegistryAccess.ecrImagePullerRole.principalArn' \
  --output text)
echo "Lightsail role: ${LIGHTSAIL_ROLE}"

# IMPORTANT: Grant Lightsail permission to pull from your ECR repo.
# Without this policy, deployments will silently fail (no app logs, just "Canceled").
aws ecr set-repository-policy \
  --repository-name islets-shiny \
  --region us-east-1 \
  --policy-text "{
  \"Version\": \"2012-10-17\",
  \"Statement\": [
    {
      \"Sid\": \"AllowLightsailPull\",
      \"Effect\": \"Allow\",
      \"Principal\": {
        \"AWS\": \"${LIGHTSAIL_ROLE}\"
      },
      \"Action\": [
        \"ecr:BatchGetImage\",
        \"ecr:GetDownloadUrlForLayer\"
      ]
    }
  ]
}"
```

Replace <YOUR_ACCOUNT_ID> with actual AWS account ID from the first commando rcapture in bash as shown above.


## Create Deployment Configuration

Create `containers.json` (update the image version number) to use the ECR image URL:

```json
{
  "islets-shiny": {
    "image": "<YOUR_ACCOUNT_ID>.dkr.ecr.us-east-1.amazonaws.com/islets-shiny:latest",
    "ports": {
      "3838": "HTTP"
    }
  }
}
```

Create public-endpoint.json:


```json
{
    "containerName": "islets-shiny",
    "containerPort": 3838,
    "healthCheck": {
        "path": "/",
        "intervalSeconds": 60,
        "timeoutSeconds": 10,
        "healthyThreshold": 2,
        "unhealthyThreshold": 5,
        "successCodes": "200-499"
    }
}

}
```
## Deploy the Application

```sh
aws lightsail create-container-service-deployment \
  --service-name islets-app \
  --containers file://containers.json \
  --public-endpoint file://public-endpoint.json
```

## Monitor Deployment

```sh
# Check deployment status
aws lightsail get-container-services --service-name islets-app

# Get the public URL once deployed
aws lightsail get-container-services \
  --service-name islets-app \
  --query 'containerServices[0].url' \
  --output text
```


Your app will be available at: https://islets-app.[random].[region].cs.amazonlightsail.com/

### Updating the Application

To deploy updates:

```sh
# 1. Rebuild the image
docker build -t islets-shiny .

# 2. Push new version
AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)
echo ${AWS_ACCOUNT_ID}

# Create ECR repository
aws ecr create-repository --repository-name islets-shiny --region us-east-1

# Login to ECR
aws ecr get-login-password --region us-east-1 | docker login --username AWS --password-stdin ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com

# Tag your image for ECR
docker tag islets-shiny:latest ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/islets-shiny:latest

# Push to ECR (this may take 10-20 minutes for a 2.5GB image)
docker push ${AWS_ACCOUNT_ID}.dkr.ecr.us-east-1.amazonaws.com/islets-shiny:latest

# Enable private registry access for Lightsail
aws lightsail update-container-service \
  --service-name islets-app \
  --private-registry-access '{"ecrImagePullerRole":{"isActive":true}}'

# 3. Update containers.json with new version number (e.g., .2, .3, etc.)

# 4. Redeploy
aws lightsail create-container-service-deployment \
  --service-name islets-app \
  --containers file://containers.json \
  --public-endpoint file://public-endpoint.json
```


### Scaling Options

If the app needs more resources:

```sh
# Upgrade to large (4GB RAM, 2 vCPU) - $80/month
aws lightsail update-container-service \
  --service-name islets-app \
  --power large

# Available power options:
# - micro: 512MB, 0.25 vCPU ($10/month)
# - small: 1GB, 0.5 vCPU ($25/month)
# - medium: 2GB, 1 vCPU ($40/month) - RECOMMENDED START
# - large: 4GB, 2 vCPU ($80/month)
# - xlarge: 8GB, 4 vCPU ($160/month)
```

### Monitoring

```sh
# View container logs
aws lightsail get-container-log \
  --service-name islets-app \
  --container-name islets-shiny

# View service metrics in AWS Console
# https://lightsail.aws.amazon.com/
```

### Check deployment status and logs

```bash
# Quick status check (state, power, deployment status)
aws lightsail get-container-services --service-name islets-app \
  --query 'containerServices[0].{state:state,power:power,currentDeployment:currentDeployment.state,nextDeployment:nextDeployment.state}' \
  --output json

# View container logs
aws lightsail get-container-log \
  --service-name islets-app \
  --container-name islets-shiny

# View deployment history (shows all past deployments and their states)
aws lightsail get-container-service-deployments \
  --service-name islets-app

# Full service details
aws lightsail get-container-services --service-name islets-app
```

### Deployments fail silently (no app logs, just "Started node" then "Canceled")

This means Lightsail can't pull the image from ECR. Two things must both be set:

1. Lightsail's ECR puller role must be active:
```bash
aws lightsail update-container-service \
  --service-name islets-app \
  --private-registry-access '{"ecrImagePullerRole":{"isActive":true}}'
```

2. Your ECR repo must have a policy granting that role pull access:
```bash
LIGHTSAIL_ROLE=$(aws lightsail get-container-services \
  --service-name islets-app \
  --query 'containerServices[0].privateRegistryAccess.ecrImagePullerRole.principalArn' \
  --output text)

aws ecr set-repository-policy \
  --repository-name islets-shiny \
  --region us-east-1 \
  --policy-text "{
  \"Version\": \"2012-10-17\",
  \"Statement\": [
    {
      \"Sid\": \"AllowLightsailPull\",
      \"Effect\": \"Allow\",
      \"Principal\": { \"AWS\": \"${LIGHTSAIL_ROLE}\" },
      \"Action\": [\"ecr:BatchGetImage\", \"ecr:GetDownloadUrlForLayer\"]
    }
  ]
}"
```

Then redeploy:
```bash
aws lightsail create-container-service-deployment \
  --service-name islets-app \
  --containers file://containers.json \
  --public-endpoint file://public-endpoint.json
```

### App won't start (logs show R errors)
* Verify DATA/Islets4-slimmed.Rds exists in the Docker image
* Ensure all R packages installed successfully

### Out of memory errors
* `medium` (2GB) works for this app, but upgrade to `large` (4GB) if you see OOM crashes
* Consider further reducing the Seurat object size

### Push command hangs
* Image is large (2.5GB), allow 15-20 minutes
* Check network connection
* Try --debug flag for verbose output
