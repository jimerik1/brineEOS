# config/__init__.py
"""
Import all salt configurations for easy access
"""
from config.nacl_config import CONFIG as NACL_CONFIG
from config.kcl_config import CONFIG as KCL_CONFIG
from config.cacl2_config import CONFIG as CACL2_CONFIG
from config.cabr2_config import CONFIG as CABR2_CONFIG
from config.znbr2_config import CONFIG as ZNBR2_CONFIG
from config.zncl2_config import CONFIG as ZNCL2_CONFIG

# Dictionary of all salt configurations
SALT_CONFIGS = {
    'NaCl': NACL_CONFIG,
    'KCl': KCL_CONFIG,
    'CaCl2': CACL2_CONFIG,
    'CaBr2': CABR2_CONFIG,
    'ZnBr2': ZNBR2_CONFIG,
    'ZnCl2': ZNCL2_CONFIG
}