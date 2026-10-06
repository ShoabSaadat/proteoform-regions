"""Smoke test: package imports and exposes its version."""


def test_import_and_version():
    import proteoform_regions

    assert proteoform_regions.__version__ == "0.1.0"


def test_module_surface():
    import proteoform_regions as pfr

    # P2 will populate these; the version check is the only hard assertion in P1.
    assert hasattr(pfr, "__version__")
