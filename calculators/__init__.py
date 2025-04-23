# calculators/__init__.py
"""
Factory function to get the appropriate calculator
"""
from calculators.single_brine_calculator import SingleBrineCalculator
from calculators.mixed_brine_calculator import MixedBrineCalculator
from config import SALT_CONFIGS

def get_calculator(brine_type, salt_composition=None):
    """
    Get the appropriate calculator for the brine type.
    
    Args:
        brine_type (str): Type of brine ('single', 'mixed', or a specific salt name)
        salt_composition (dict, optional): Composition of mixed brine
        
    Returns:
        Calculator: Instance of the appropriate calculator
        
    Raises:
        ValueError: If brine_type is not supported
    """
    if brine_type.lower() == 'mixed':
        return MixedBrineCalculator(SALT_CONFIGS)
    
    elif brine_type.lower() == 'single':
        # For 'single' type, we need to determine which salt from the composition
        if not salt_composition or len(salt_composition) != 1:
            raise ValueError("For 'single' brine type, exactly one salt must be specified in salt_composition")
        
        salt_name = next(iter(salt_composition.keys()))
        if salt_name not in SALT_CONFIGS:
            raise ValueError(f"Unsupported salt: {salt_name}. Supported salts are: {', '.join(SALT_CONFIGS.keys())}")
        
        return SingleBrineCalculator(SALT_CONFIGS[salt_name])
    
    else:
        # Check if brine_type is a specific salt name
        if brine_type in SALT_CONFIGS:
            return SingleBrineCalculator(SALT_CONFIGS[brine_type])
        
        # If we get here, the brine type is not supported
        raise ValueError(f"Unsupported brine type: {brine_type}. "
                         f"Supported types are: 'single', 'mixed', or one of {', '.join(SALT_CONFIGS.keys())}")