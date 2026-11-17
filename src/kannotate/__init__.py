"""Kannotate's file-independent protein annotation API."""

from kbiosuite.version import __version__

from .engine import AnnotationEngine
from .models import Annotation, Hit, Protein

__all__ = ["Annotation", "AnnotationEngine", "Hit", "Protein", "__version__"]
