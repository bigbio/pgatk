"""Unit tests for the COSMIC downloader URL construction and auth token.

No network. Two tests:
  * `test_build_api_url_pattern` — verifies build_api_url() composes the
    scripted-download endpoint correctly from server, api_endpoint, path
    and bucket.
  * `test_basic_auth_token_format` — verifies the base64(user:password)
    token has no trailing newline (a common pitfall when echo'ing).
"""
from pgatk.cgenomes.cosmic_downloader import CosmicDownloadService
from pgatk.config.registry import load_config


def test_build_api_url_pattern(tmp_path):
    """The API URL should be {server}/{endpoint}?path={path}&bucket={bucket}."""
    config_data = load_config("cosmic", None)
    pipeline_arguments = {
        CosmicDownloadService.CONFIG_OUTPUT_DIRECTORY: str(tmp_path),
        CosmicDownloadService.CONFIG_COSMIC_FTP_USER: "user@example.com",
        CosmicDownloadService.CONFIG_COSMIC_FTP_PASSWORD: "secret",
    }
    svc = CosmicDownloadService(config_data, pipeline_arguments)

    url = svc.build_api_url("grch38/cosmic/v103/Cosmic_MutantCensus_Tsv_v103_GRCh38.tar")

    assert url == (
        "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted"
        "?path=grch38/cosmic/v103/Cosmic_MutantCensus_Tsv_v103_GRCh38.tar"
        "&bucket=downloads"
    )


def test_basic_auth_token_format(tmp_path):
    """Auth token is base64(user:password), no trailing newline."""
    import base64

    config_data = load_config("cosmic", None)
    pipeline_arguments = {
        CosmicDownloadService.CONFIG_OUTPUT_DIRECTORY: str(tmp_path),
        CosmicDownloadService.CONFIG_COSMIC_FTP_USER: "user@example.com",
        CosmicDownloadService.CONFIG_COSMIC_FTP_PASSWORD: "secret",
    }
    svc = CosmicDownloadService(config_data, pipeline_arguments)

    expected = base64.b64encode(b"user@example.com:secret").decode("utf-8")
    assert svc._cosmic_token == expected
