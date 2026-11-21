"""HTTPS trust for the suite's stdlib urllib clients: Kmatch's RCSB access, Kannotate's downloads.

macOS decides which certificates to trust. `truststore` hands each server's chain to
Security.framework, which evaluates it for the requested host name with the trust settings
that apply to this application: Apple's current system roots, plus roots that an administrator
or the user has trusted. The decision follows Apple's trust updates.

Without this, the Python interpreter the apps ship would use OpenSSL's default file,
`/private/etc/ssl/cert.pem`. That is a separate list, not kept in step with the system roots:
a public host that moves to a newer root fails there while macOS itself still trusts it.

There is deliberately no fallback to that file or to a bundled CA list, and `SSL_CERT_FILE`
and `SSL_CERT_DIR` do not change the decision. Kblast's `requests` client is separate.
"""

import ssl
import urllib.request


def context():
    """A new client context; certificate and hostname verification are always on."""
    import truststore  # first HTTPS use, never application launch

    made = truststore.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
    made.set_alpn_protocols(["http/1.1"])  # http.client only does this for its own default
    if made.verify_mode != ssl.CERT_REQUIRED or not made.check_hostname:
        raise RuntimeError("TLS certificate and hostname verification must be enabled")
    return made


class NoRedirect(urllib.request.HTTPRedirectHandler):
    """A redirect becomes its own HTTP error instead of a request to another address."""

    def redirect_request(self, *arguments, **options):
        return None


def urlopen(request, *, timeout, direct=False):
    """`urllib.request.urlopen` with macOS trust.

    direct=True is for qualification. It disables proxies and redirects.
    HTTPS requests use the verified context for the requested host.
    """
    # One context per call: truststore serialises handshakes that share a context.
    handlers = [urllib.request.HTTPSHandler(context=context())]
    if direct:
        handlers += [urllib.request.ProxyHandler({}), NoRedirect()]
    return urllib.request.build_opener(*handlers).open(request, timeout=timeout)


def head(url):
    """Explicit qualification only: one direct HEAD request, returning its HTTP status."""
    request = urllib.request.Request(url, method="HEAD")
    with urlopen(request, timeout=20, direct=True) as response:
        return response.status
