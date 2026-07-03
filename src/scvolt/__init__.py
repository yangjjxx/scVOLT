"""Public package metadata for scVOLT."""

from importlib import import_module

from scvolt.model import SCVOLTModel

__all__ = ["SCVOLTModel", "tools"]
__version__ = "1.0.0"


def __getattr__(name):
    """Lazily expose public package objects listed in ``__all__``."""
    if name == "tools":
        value = import_module(".tools", __name__)
    elif name == "scVOLTModel":
        # Backward-compatible alias for older notebooks that used this spelling.
        value = SCVOLTModel
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    globals()[name] = value
    return value
