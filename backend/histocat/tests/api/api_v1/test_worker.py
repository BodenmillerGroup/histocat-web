import requests

from histocat.config import config
from histocat.tests.utils.utils import get_server_api


def test_worker_test(superuser_token_headers):
    server_api = get_server_api()
    data = {"msg": "test"}
    r = requests.post(
        f"{server_api}{config.API_V1_STR}/utils/test-worker/",
        json=data,
        headers=superuser_token_headers,
    )
    response = r.json()
    assert response["msg"] == "Word received"
