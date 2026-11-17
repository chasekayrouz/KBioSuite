"""One local operation owns cancellation, publication and rollback of new cache data."""

import threading
from contextlib import contextmanager
from contextvars import ContextVar

_current = ContextVar("kannotate_run", default=None)


class RunCancelled(Exception):
    """An intentional stop, never a failed search or an empty scientific result."""


class RunControl:
    def __init__(self):
        self._stopped = threading.Event()
        self._lock = threading.RLock()
        self._finished = False
        self._rollback = []

    def request_stop(self, before_stop=None):
        with self._lock:
            if self._finished:
                return False
            if before_stop:
                before_stop()
            self._stopped.set()
            return True

    def check(self):
        if self._stopped.is_set():
            raise RunCancelled("Search stopped")

    @contextmanager
    def publication(self):
        # A completed output wins this lock before Stop, or Stop prevents publication.
        with self._lock:
            self.check()
            yield

    def on_rollback(self, callback):
        with self._lock:
            self._rollback.append(callback)

    def finish(self, publish=None):
        with self.publication():
            if publish:
                publish()
            self._finished = True
            self._rollback.clear()

    def rollback(self):
        errors = []
        for callback in reversed(self._rollback):
            try:
                callback()
            except OSError as exc:
                errors.append(str(exc))
        self._rollback.clear()
        if errors:
            raise OSError("Could not remove all interrupted search data: " + "; ".join(errors))

    @contextmanager
    def activate(self):
        """Make this control available to adapters, including in engine worker threads."""
        token = _current.set(self)
        try:
            self.check()
            yield self
        finally:
            _current.reset(token)


def current_control():
    return _current.get()


def check_cancelled():
    control = current_control()
    if control is not None:
        control.check()


@contextmanager
def run_scope(control=None):
    """Only the outer operation commits/rolls back; nested APIs share its control."""
    inherited = current_control()
    control = control or inherited or RunControl()
    if control is inherited:
        control.check()
        yield control
        return
    with control.activate():
        try:
            yield control
            control.finish()
        except BaseException:
            control.rollback()
            raise


class CheckedText:
    """Check between format-library reads/writes without changing serialized text."""

    def __init__(self, handle):
        self.handle = handle

    def __getattr__(self, name):
        return getattr(self.handle, name)

    def read(self, *args):
        check_cancelled()
        return self.handle.read(*args)

    def readline(self, *args):
        check_cancelled()
        return self.handle.readline(*args)

    def write(self, text):
        check_cancelled()
        return self.handle.write(text)

    def __iter__(self):
        return self

    def __next__(self):
        check_cancelled()
        return next(self.handle)
