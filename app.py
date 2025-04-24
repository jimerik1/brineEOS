# app.py
"""
Main Flask application
"""
import os
from flask import Flask
from routes.routes import api_bp
from routes.water_routes import water_api_bp # Import the new water blueprint
import logging
from logging.handlers import RotatingFileHandler
from flask_cors import CORS

def create_app():
    """
    Create and configure the Flask application.

    Returns:
        app: Configured Flask application
    """
    app = Flask(__name__)

    # Enable Cross-Origin Resource Sharing
    CORS(app)

    # Configure logging
    if not app.debug:
        if not os.path.exists('logs'):
            os.mkdir('logs')
        file_handler = RotatingFileHandler('logs/brine_density_api.log', maxBytes=10240, backupCount=10)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'
        ))
        file_handler.setLevel(logging.INFO)
        app.logger.addHandler(file_handler)
        app.logger.setLevel(logging.INFO)
        app.logger.info('Brine Density API startup')

    # Register blueprints
    app.register_blueprint(api_bp)
    app.register_blueprint(water_api_bp) # Register the new water API blueprint

    # Health check endpoint
    @app.route('/health', methods=['GET'])
    def health_check():
        return {'status': 'ok'}

    # Error handlers
    @app.errorhandler(404)
    def not_found_error(error):
        return {'error': 'Resource not found'}, 404

    @app.errorhandler(500)
    def internal_error(error):
        app.logger.error('Server Error: %s', error)
        return {'error': 'Internal server error'}, 500

    return app

if __name__ == "__main__":
    app = create_app()
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)
