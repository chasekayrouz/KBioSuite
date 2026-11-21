import sys

if __name__ == "__main__":
    from kbiosuite.runtime_checks import dispatch_network_check

    def probes():
        from kannotate.databases import network_probes

        return network_probes()

    if (status := dispatch_network_check("Kannotate", "urllib", probes)) is not None:
        raise SystemExit(status)
    if len(sys.argv) == 3 and sys.argv[1] == "--kannotate-resident-worker":
        from kannotate.search.resident import worker_main

        worker_main(sys.argv[2])
    else:
        from kannotate.gui import main

        main()
