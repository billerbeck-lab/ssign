# Vendored from https://github.com/zhengdd0422/PLM-Effector's upstream
# `trainer.py` (CC-BY 3.0). Renamed to `inference.py` because the kept
# functions are all inference-time — there is no training code here.
# Upstream imports of `random`, `numpy`, and the large block of
# sklearn.metrics have been dropped; only `os`/`gc`/`torch` remain since
# the prediction helpers below memory-map model files (os.path.join) and
# explicitly free GPU memory between models (gc.collect).
import gc
import os

import torch
from torch.utils.data import DataLoader, TensorDataset


def test_4predict_inbatch(model, features, device, batch_size=32):
    model.eval()
    dataset = TensorDataset(features)
    loader = DataLoader(dataset, batch_size=batch_size)

    all_preds = []
    all_probs = []

    with torch.no_grad():
        for batch in loader:
            x = batch[0].to(device)
            outputs = model(x).squeeze(dim=1)
            probs = torch.sigmoid(outputs)
            preds = (probs > 0.5).long()

            all_preds.append(preds.cpu())
            all_probs.append(probs.cpu())

    all_preds = torch.cat(all_preds).numpy()
    all_probs = torch.cat(all_probs).numpy()
    return all_preds, all_probs


def loadmodel_4predict(model_dir, model_name, x_test, device):
    model = torch.load(os.path.join(model_dir, model_name), map_location=device)
    if isinstance(model, torch.nn.DataParallel):
        model = model.module
        model.eval()
    _, test_probs = test_4predict_inbatch(model, x_test, device, batch_size=128)
    del model
    torch.cuda.empty_cache()
    gc.collect()
    return test_probs


def loadmodel_4test(model_dir, model_name, x_test, device):
    model = torch.load(os.path.join(model_dir, model_name), map_location=device)
    if isinstance(model, torch.nn.DataParallel):
        model = model.module
        model.eval()
    test_preds, test_probs = test_4predict_inbatch(
        model, x_test, device, batch_size=128
    )
    del model
    torch.cuda.empty_cache()
    gc.collect()

    return test_preds, test_probs
