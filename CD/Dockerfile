# Dockerfile
FROM python:3.9-slim

WORKDIR /app

# Copy requirements and install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# add src to path
ENV PYTHONPATH=/app/src
ENV PYTHONPATH=$PYTHONPATH

# Environment variables
ENV FLASK_APP=src/app.py
ENV FLASK_ENV=production
ENV PYTHONUNBUFFERED=1

# Expose port
EXPOSE 5000

# Run gunicorn server
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "--workers", "4", "app:create_app()"]