
"""
Toolbox related to REST services
"""

import requests
from ratelimit import limits, sleep_and_retry

ONE_MINUTE = 60

@sleep_and_retry
@limits(calls=1000, period=ONE_MINUTE)
def call_api(url):
    response = requests.get(url, headers={"Content-Type": 'application/json'})

    if response.status_code != 200:
        raise Exception('API response: {}'.format(response.status_code))
    return response

@sleep_and_retry
@limits(calls=1000, period=ONE_MINUTE)
def call_api_raw(url):
    response = requests.get(url)

    if response.status_code != 200:
        raise Exception('API response: {}'.format(response.status_code))
    return response
def make_rest_request(url):
    response = requests.get(url, headers={"Content-Type": "application/json"})
    if not response.ok:
        response.raise_for_status()
    return response.json()


if __name__ == '__main__':
    print("ERROR: This script is part of a pipeline collection and it is not meant to be run in stand alone mode")
