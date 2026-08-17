import requests

# GitHub Search API endpoint
url = "https://api.github.com/search/code"
params = {
    "q": "extension:ipynb",
    "per_page": 1
}

# Optional: Add your GitHub personal access token to increase rate limits
headers = {
    "Accept": "application/vnd.github.v3+json"
}

response = requests.get(url, params=params, headers=headers)
data = response.json()

if "total_count" in data:
    print(f"Approximate number of Jupyter notebooks on GitHub: {data['total_count']:,}")
else:
    print("Could not retrieve count:", data)
