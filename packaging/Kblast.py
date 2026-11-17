from kblast.app import main

if __name__ == "__main__":
    # The only entry point that imports the requests probe, so only Kblast bundles requests.
    from kbiosuite.requests_probe import head
    from kbiosuite.runtime_checks import dispatch_network_check

    if (status := dispatch_network_check("Kblast", "requests", head=head)) is not None:
        raise SystemExit(status)
    main()
