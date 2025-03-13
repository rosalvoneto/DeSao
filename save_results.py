import json
import os

def save_results(dic_temp, dic_delta, dic_prob, dic_accept, dic_news_ok, dic_news, results_dir="Results"):
    """
    Saves Simulated Annealing results to JSON files.
    
    Args:
        dic_temp (dict): Dictionary containing temperatures
        dic_delta (dict): Dictionary containing deltas
        dic_prob (dict): Dictionary containing probabilities
        dic_accept (dict): Dictionary containing acceptances
        dic_news_ok (dict): Dictionary containing valid new molecules
        dic_news (dict): Dictionary containing all new molecules
        results_dir (str): Directory where results will be saved
    """
    # Ensure directory exists
    os.makedirs(results_dir, exist_ok=True)
    
    results = {
        "temperatures.json": dic_temp,
        "deltas.json": dic_delta,
        "probabilities.json": dic_prob,
        "accepted.json": dic_accept,
        "dic_news_ok.json": dic_news_ok,
        "dic_news.json": dic_news
    }
    
    for filename, data in results.items():
        filepath = os.path.join(results_dir, filename)
        with open(filepath, 'w') as file:
            json.dump(data, file) 