import re
from pathlib import Path
from unittest.mock import Mock

import pytest
import webview

import kannotate
from kannotate.gui import AnnotationJob, DesktopApi, create_app


def desktop_api(result):
    api = DesktopApi(AnnotationJob())
    api._window = Mock()
    api._window.create_file_dialog.return_value = result
    return api


@pytest.mark.parametrize("result", [None, (), ("/tmp/input with spaces.faa", "/tmp/cluster.gbk")])
def test_input_picker_uses_pinned_enum_multi_select_and_filters(result):
    api = desktop_api(result)
    assert api.choose_inputs() == list(result or [])
    args, kwargs = api._window.create_file_dialog.call_args
    assert args == (webview.FileDialog.OPEN,)
    assert kwargs["allow_multiple"] is True
    assert "*.faa" in kwargs["file_types"][0] and "*.gbk" in kwargs["file_types"][0]
    assert not hasattr(api, "window") and not hasattr(api, "job")


@pytest.mark.parametrize("result", [None, (), ("/tmp/output with spaces",)])
def test_output_picker_handles_cancel_and_folder_result(result):
    api = desktop_api(result)
    assert api.choose_output() == (result[0] if result else None)
    api._window.create_file_dialog.assert_called_once_with(webview.FileDialog.FOLDER)


def test_missing_window_and_native_failure_are_actionable():
    api = DesktopApi(AnnotationJob())
    with pytest.raises(RuntimeError, match="Close and reopen"):
        api.choose_inputs()
    api._window = Mock()
    api._window.create_file_dialog.side_effect = RuntimeError("Native picker unavailable")
    with pytest.raises(RuntimeError, match="Native picker"):
        api.choose_output()


def test_reveal_only_exposes_current_successful_outputs(monkeypatch):
    api = desktop_api(None)
    call = Mock()
    monkeypatch.setattr("kannotate.gui.subprocess.run", call)
    with pytest.raises(ValueError, match="output from this run"):
        api.reveal("/tmp/unrelated")
    call.assert_not_called()
    api._job.update(results=[{"output": "/tmp/annotation output.xlsx"}])
    api.reveal("/tmp/annotation output.xlsx")
    call.assert_called_once_with(["/usr/bin/open", "-R", "/tmp/annotation output.xlsx"], check=True)
    api._job.update(
        results=[
            {
                "output": "/tmp/annotation.gbk",
                "outputs": ["/tmp/annotation.gbk", "/tmp/annotation.xlsx"],
            }
        ]
    )
    for path in ("/tmp/annotation.gbk", "/tmp/annotation.xlsx"):
        api.reveal(path)
        call.assert_called_with(["/usr/bin/open", "-R", path], check=True)


def test_bridge_csp_is_native_only_and_buttons_start_disabled():
    for desktop in [False, True]:
        response = create_app(AnnotationJob(), desktop=desktop).test_client().get("/")
        policy = response.headers["Content-Security-Policy"]
        assert ("'unsafe-eval'" in policy) is desktop
        assert "'unsafe-inline'" not in policy
        assert b'id="choose-inputs" disabled' in response.data
        assert b'id="choose-output" disabled' in response.data
        assert b"/static/bridge.js" in response.data


def test_normal_close_releases_session_but_busy_close_waits():
    api = DesktopApi(AnnotationJob())
    close = Mock()
    api._job.search_session.close = close
    api._job.update(busy=True)
    assert api.can_close() is False
    close.assert_not_called()
    api._job.update(busy=False)
    assert api.can_close() is True
    close.assert_called_once()


def test_page_script_only_references_elements_in_the_page():
    """A missing id would make a button handler throw (e.g. the setup result line)."""
    static = Path(kannotate.__file__).with_name("static")
    ids = set(re.findall(r'id="([^"]+)"', (static / "index.html").read_text()))
    references = set(re.findall(r"\$\('([^']+)'\)", (static / "app.js").read_text()))
    assert "setup-error" in references and references <= ids
