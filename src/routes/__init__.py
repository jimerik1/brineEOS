# routes/__init__.py
from .routes import api_bp
from .water_routes import water_api_bp

__all__ = ["api_bp", "water_api_bp"]