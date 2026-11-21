from kmatch.gui import main

if __name__ == "__main__":
    from kbiosuite.runtime_checks import dispatch_network_check
    from kmatch.pdb_data import network_probes
    from kmatch.setup import dispatch_uv_check

    # Kmatch's own HTTPS client (RCSB access) is stdlib urllib with macOS trust.
    if (status := dispatch_network_check("Kmatch", "urllib", network_probes)) is not None:
        raise SystemExit(status)
    if (status := dispatch_uv_check()) is not None:
        raise SystemExit(status)
    main()
