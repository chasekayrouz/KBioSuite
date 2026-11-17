from kmatch.gui import main

if __name__ == "__main__":
    from kbiosuite.runtime_checks import dispatch_network_check
    from kmatch.setup import dispatch_uv_check

    # Kmatch's own HTTPS client (RCSB access) is stdlib urllib, like Kannotate's.
    if (status := dispatch_network_check("Kmatch", "urllib")) is not None:
        raise SystemExit(status)
    if (status := dispatch_uv_check()) is not None:
        raise SystemExit(status)
    main()
