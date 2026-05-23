{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "33b087c4-5902-4d34-b649-fa20a4660689",
   "metadata": {},
   "outputs": [],
   "source": [
    "#import\n",
    "import pandas as pd\n",
    "from rdkit import Chem\n",
    "from rdkit.Chem import Descriptors\n",
    "from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator\n",
    "import numpy as np\n",
    "\n",
    "c_ter=pd.read_csv(\"c-terminus-smiles.csv\")\n",
    "n_ter=pd.read_csv(\"n-terminus-smiles.csv\")\n",
    "uncharged_ter=pd.read_csv(\"uncharged-terminus-smiles.csv\")\n",
    "\n",
    "\n",
    "# 1) Build a calculator for ALL RDKit 1D+2D descriptors\n",
    "desc_names = [name for name, _ in Descriptors.descList]   # typically ~200 2D descriptors\n",
    "calc = MolecularDescriptorCalculator(desc_names)\n",
    "\n",
    "def smiles_to_descriptors(smiles: str):\n",
    "    \"\"\"Return a dict {descriptor_name: value} for one SMILES.\n",
    "       Invalid/empty SMILES -> all NaN.\"\"\"\n",
    "    if not isinstance(smiles, str) or not smiles.strip():\n",
    "        return {n: np.nan for n in desc_names}\n",
    "    mol = Chem.MolFromSmiles(smiles)\n",
    "    if mol is None:\n",
    "        return {n: np.nan for n in desc_names}\n",
    "    try:\n",
    "        vals = calc.CalcDescriptors(mol)\n",
    "        return dict(zip(desc_names, vals))\n",
    "    except Exception:\n",
    "        return {n: np.nan for n in desc_names}\n",
    "\n",
    "# 2) Compute descriptor table\n",
    "c_des = c_ter['SMILES'].apply(smiles_to_descriptors).apply(pd.Series)\n",
    "c_des=c_des.drop(\"Ipc\",axis=1)\n",
    "n_des = n_ter['SMILES'].apply(smiles_to_descriptors).apply(pd.Series)\n",
    "n_des=n_des.drop(\"Ipc\",axis=1)\n",
    "uncharged_des = uncharged_ter['SMILES'].apply(smiles_to_descriptors).apply(pd.Series)\n",
    "uncharged_des=uncharged_des.drop(\"Ipc\",axis=1)\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "\n",
    "# 3) Combine with your original data\n",
    "c_des = pd.concat([c_ter.reset_index(drop=True), c_des], axis=1)\n",
    "n_des = pd.concat([n_ter.reset_index(drop=True), n_des], axis=1)\n",
    "uncharged_des = pd.concat([uncharged_ter.reset_index(drop=True), uncharged_des], axis=1)\n",
    "\n",
    "\n",
    "\n",
    "# 4) (optional) save\n",
    "c_des.to_csv(\"c-terminus-descriptors.csv\", index=False)\n",
    "n_des.to_csv(\"n-terminus-descriptors.csv\", index=False)\n",
    "uncharged_des.to_csv(\"uncharged-terminus-descriptors.csv\", index=False)\n",
    "\n"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.10"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
