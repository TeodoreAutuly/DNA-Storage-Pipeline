# DNA Storage Pipeline

Un pipeline complet de simulation pour le stockage de données sur ADN synthétique, incluant encodage, simulation de canal bruité, décodage par consensus et correction d'erreurs LDPC.

---

## Table des matières

1. [Vue d'ensemble](#vue-densemble)
2. [Fonctionnalités](#fonctionnalités)
3. [Architecture du projet](#architecture-du-projet)
4. [Installation](#installation)
5. [Utilisation](#utilisation)
6. [Modèles statistiques](#modèles-statistiques)
7. [Interface web interactive](#interface-web-interactive)
8. [Technologies utilisées](#technologies-utilisées)
9. [Contexte scientifique](#contexte-scientifique)
10. [Licence](#licence)
11. [Auteur](#auteur)

---

## Vue d'ensemble

Ce projet implémente un **simulateur de stockage de données sur ADN**, permettant de transformer des données numériques en séquences d'ADN synthétique, de simuler les erreurs du processus physique (synthèse, stockage, séquençage), puis de reconstruire les données originales grâce à des algorithmes de correction d'erreurs.

### Principe général

```
Données binaires
    ↓
Encodage LDPC (ajout de redondance)
    ↓
Conversion en séquence ADN (A, C, G, T)
    ↓
Simulation du canal bruité (erreurs de synthèse/séquençage)
    ↓
Algorithme de consensus (alignement de lectures multiples)
    ↓
Décodage LDPC (correction d'erreurs)
    ↓
Données binaires restaurées
```

### Objectifs du projet

- **Simuler fidèlement** les erreurs réelles observées lors de la synthèse et du séquençage d'ADN
- **Évaluer la robustesse** des codes correcteurs d'erreurs (LDPC)
- **Fournir une interface pédagogique** pour comprendre le stockage ADN
- **Supporter différents modèles** de canaux de séquençage (k4_idt_0124, k6_dnarxiv)

---

## Fonctionnalités

### Simulation complète du pipeline

- **Encodage binaire vers ADN** : conversion des données en séquences quaternaires (A, C, G, T)
- **Code LDPC paramétrable** : protection contre les erreurs avec ratio d'information configurable
- **Modèles de canal réalistes** : simulation d'erreurs de substitution, insertion et délétion
- **Algorithme de consensus** : reconstruction à partir de lectures multiples et bruitées
- **Correction d'erreurs** : décodage LDPC pour restaurer les données originales

### Modèles statistiques avancés

- **Erreurs dépendantes de la position** : début, milieu et fin de séquence
- **Erreurs dépendantes du contexte** : probabilités conditionnées par les k-mers
- **Profils de qualité réalistes** : basés sur des données expérimentales (IDT, DNArchive)
- **Support multi-modèles** : k4_idt_0124 et k6_dnarxiv inclus

### Interface utilisateur

- **Interface web interactive** : visualisation en temps réel du pipeline complet
- **Paramètres configurables** : taux d'erreur, nombre de lectures, ratio LDPC
- **Visualisation pas à pas** : affichage détaillé de chaque étape du processus
- **Résultats comparatifs** : analyse du message original vs décodé

---

## Architecture du projet

```
DNA-Storage-Pipeline/
│
├── index.html                      # Interface web interactive
│
├── Models/                         # Modèles statistiques d'erreurs
│   ├── k4_idt_0124/               # Modèle IDT avec k-mers de taille 4
│   │   ├── Kval.txt               # Valeur de k
│   │   ├── A2C_Avg.txt            # Probabilités de substitution A→C
│   │   ├── AErrAvg.txt            # Taux d'erreur pour base A
│   │   ├── BegDelByPosAvg.txt     # Délétions en début de séquence
│   │   └── ...                    # Autres matrices de probabilités
│   │
│   └── k6_dnarxiv/                # Modèle DNArchive avec k-mers de taille 6
│       └── ...                     # Mêmes types de fichiers
│
├── Simulator_Python/               # Implémentation Python du simulateur
│   ├── test_simulator.ipynb       # Notebook de test et démonstration
│   ├── fastaSequences/             # Séquences FASTA pour les tests
│   └── __pycache__/               # Cache Python
│
├── src/                            # Code source Python
│   ├── lib_model.py               # Chargement et gestion des modèles
│   └── lib_utils.py               # Fonctions utilitaires (conversion, etc.)
│
├── requirements.txt                # Dépendances Python
├── python3.11_requirements.txt     # Instructions d'installation Python 3.11
├── install_muscle.txt              # Installation de MUSCLE (alignement)
│
└── LICENSE                         # Licence MIT
```

### Fichiers clés

#### Modèles statistiques

Les dossiers `Models/` contiennent des matrices de probabilités pour :

- **Substitutions** : A2C, A2G, A2T, C2A, etc. (12 fichiers)
- **Insertions/Délétions** : AInsAvg, ADelAvg, etc. (8 fichiers)
- **Erreurs par position** : BegErrByPosAvg, MidErrByPosAvg, EndErrByPosAvg
- **Contexte k-mer** : KmerTrans_RatesAvg, KmerYi_prevYiM_RatesAvg, etc.
- **Distributions de longueur** : insLenBegRates, unmapLenEndRates

#### Code source

- **lib_model.py** : classe `Str_Model` pour charger les modèles statistiques
- **lib_utils.py** : fonctions de conversion (Dec2Dna, Dna2Dec, etc.)
- **test_simulator.ipynb** : notebook Jupyter pour tester le pipeline

---

## Installation

### Prérequis

- **Python 3.11** recommandé (compatible avec 3.8+)
- **pip** pour la gestion des paquets
- **MUSCLE** (optionnel) pour l'alignement de séquences

### Installation de Python 3.11 (Ubuntu/Debian)

```bash
sudo add-apt-repository ppa:deadsnakes/ppa -y
sudo apt update && sudo apt install python3.11 python3.11-venv python3.11-dev -y
python3.11 --version
sudo apt install python3-pip
```

### Installation des dépendances Python

```bash
python3.11 -m pip install -r requirements.txt
```

Les dépendances incluent :

- **numpy** ≥ 1.20.0 : calculs numériques et matrices
- **matplotlib** ≥ 3.5.0 : visualisation des résultats
- **Pillow** ≥ 9.0.0 : traitement d'images
- **biopython** ≥ 1.79 : manipulation de séquences biologiques

### Installation de MUSCLE (optionnel)

Pour l'alignement multiple de séquences :

```bash
sudo apt install muscle
```

---

## Utilisation

### Interface web interactive

L'interface web `index.html` permet de tester le pipeline sans installation :

1. Ouvrir `index.html` dans un navigateur moderne
2. Entrer un message texte à encoder
3. Configurer les paramètres :
   - **Ratio LDPC** : 0.5 à 0.9 (plus élevé = plus dense, moins robuste)
   - **Nombre de lectures** : 5 à 50 (plus élevé = meilleur consensus)
   - **Taux d'erreur** : 0% à ~13% (simulé par le slider)
4. Cliquer sur "Lancer la Simulation"
5. Observer les 8 étapes du pipeline en temps réel

### Utilisation Python

#### Notebook Jupyter

```bash
cd Simulator_Python
jupyter notebook test_simulator.ipynb
```

#### Script Python personnalisé

```python
import sys
sys.path.append('./src')
from lib_model import Str_Model
from lib_utils import Dec2Dna

# Charger un modèle
model = Str_Model('./Models/k4_idt_0124')

# Convertir des données en ADN
data = [0, 1, 2, 3, 0, 1]  # Représentation quaternaire
dna_sequence = Dec2Dna(data)
print(dna_sequence)  # ['A', 'C', 'G', 'T', 'A', 'C']
```

---

## Modèles statistiques

### k4_idt_0124

Modèle basé sur les données de **Integrated DNA Technologies (IDT)**, janvier 2024 :

- Longueur de k-mer : **4**
- Technologie : Synthèse par phosphoramidite
- Séquençage : Illumina
- Caractéristiques : Erreurs modérées, bon pour prototypage

### k6_dnarxiv

Modèle basé sur les données de **DNArchive** :

- Longueur de k-mer : **6**
- Technologie : Synthèse avancée
- Séquençage : Oxford Nanopore / PacBio
- Caractéristiques : Plus de contexte, profil d'erreur différent

### Types d'erreurs modélisés

#### Substitutions

Remplacement d'une base par une autre (ex : A → C)

- Probabilités dépendantes de la base originale
- 12 matrices : A2C, A2G, A2T, C2A, C2G, C2T, G2A, G2C, G2T, T2A, T2C, T2G

#### Insertions

Ajout d'une ou plusieurs bases non désirées

- Probabilités par base (A, C, G, T)
- Distribution de longueur d'insertion

#### Délétions

Perte d'une ou plusieurs bases

- Probabilités par base (A, C, G, T)
- Distribution de longueur de délétion

#### Erreurs positionnelles

Les erreurs varient selon la position dans la séquence :

- **Début** (BegErrByPosAvg) : souvent plus d'erreurs
- **Milieu** (MidErrByPosAvg) : taux d'erreur stable
- **Fin** (EndErrByPosAvg) : peut augmenter légèrement

#### Contexte k-mer

Les erreurs dépendent du contexte local (k bases précédentes) :

- **KmerTrans_RatesAvg** : transitions entre k-mers
- **KmerYi_prevYiM_RatesAvg** : probabilités conditionnées par l'état précédent

---

## Interface web interactive

### Fonctionnalités de l'interface

L'interface `index.html` offre une expérience pédagogique complète :

#### Simulateur interactif

- Entrée de texte libre
- Contrôles en temps réel :
  - Curseur de ratio LDPC (0.5 - 0.9)
  - Curseur de nombre de lectures (5 - 50)
  - Curseur de taux d'erreur (0% - 13%)
- Bouton de lancement avec effet lumineux

#### Visualisation du pipeline

Affichage séquentiel des 8 étapes :

1. **Encodage binaire** : texte → bits
2. **Protection LDPC** : ajout de bits de parité
3. **Encodage ADN** : binaire → A/C/G/T
4. **Canal bruité** : simulation d'erreurs multiples
5. **Consensus** : vote majoritaire sur lectures
6. **Décodage binaire** : ADN → binaire synchronisé
7. **Correction LDPC** : utilisation des bits de parité
8. **Résultat final** : comparaison original/décodé

#### Informations éducatives

- **4 étapes physiques** : synthèse, stockage, séquençage, décodage
- **Avantages du stockage ADN** :
  - Densité extrême (200 exaoctets/gramme)
  - Durabilité millénaire
- **Défis actuels** :
  - Coût élevé
  - Vitesse d'écriture/lecture

### Design

- **Thème sombre** : confort visuel pour code et données
- **Animations fluides** : transitions CSS3
- **Responsive** : adapté mobile, tablette, desktop
- **Typographie** : famille Inter pour modernité
- **Couleurs** : palette bleu nuit/cyan pour cohérence scientifique

---

## Technologies utilisées

### Backend / Simulation

| Technologie | Version | Utilisation |
|------------|---------|-------------|
| Python | 3.11+ | Langage principal |
| NumPy | ≥ 1.20.0 | Calculs matriciels et statistiques |
| Matplotlib | ≥ 3.5.0 | Visualisation de données |
| BioPython | ≥ 1.79 | Manipulation de séquences ADN |
| Pillow | ≥ 9.0.0 | Traitement d'images |
| MUSCLE | - | Alignement multiple de séquences |

### Frontend

| Technologie | Utilisation |
|------------|-------------|
| HTML5 | Structure de l'interface |
| CSS3 / Tailwind CSS | Styling responsive et moderne |
| JavaScript (Vanilla) | Logique de simulation côté client |
| Google Fonts (Inter) | Typographie |

### Algorithmes

- **Code LDPC** : Low-Density Parity-Check pour correction d'erreurs
- **Consensus** : Vote majoritaire sur lectures multiples
- **K-mer modeling** : Modèles statistiques contextuels
- **Synchronisation** : Alignement de séquences bruitées

---

## Contexte scientifique

### Le stockage sur ADN

L'ADN est le support d'information le plus dense et durable connu dans l'univers. Un gramme d'ADN peut théoriquement stocker **plus de 200 exaoctets** de données, soit l'équivalent de millions de disques durs modernes.

### Processus physique

#### 1. Synthèse (Écriture)

Les données binaires sont converties en séquences d'ADN (A, C, G, T). Des machines spécialisées assemblent chimiquement ces brins base par base.

#### 2. Stockage (Archivage)

Les molécules d'ADN sont déshydratées et encapsulées dans des billes de silice pour une conservation millénaire sans énergie.

#### 3. Séquençage (Lecture)

Les capsules sont réhydratées et analysées par un séquenceur qui lit l'ordre des bases, générant des millions de lectures bruitées.

#### 4. Décodage (Interprétation)

Des algorithmes (consensus + LDPC) reconstruisent la séquence parfaite et la reconvertissent en données numériques.

### Défis techniques

#### Erreurs de séquençage

Les technologies actuelles introduisent des erreurs :

- **Substitutions** : ~0.1% à 1% (base changée)
- **Insertions** : ~0.01% à 0.5% (base ajoutée)
- **Délétions** : ~0.01% à 0.5% (base manquante)

#### Solutions implémentées

- **Redondance** : codes correcteurs d'erreurs (LDPC, Reed-Solomon)
- **Consensus** : séquençage multiple pour vote majoritaire
- **Modèles statistiques** : prédiction des erreurs probables

### Applications futures

- **Archivage à long terme** : données historiques, culturelles
- **Sauvegarde critique** : génomes, connaissances humaines
- **Stockage spatial** : missions longue durée, résistance aux radiations

---

## Licence

Ce projet est distribué sous licence **MIT**.

```
MIT License

Copyright (c) 2025 TeodoreAutuly

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

Voir le fichier [LICENSE](LICENSE) pour plus de détails.

---

## Auteur

**Téodore Autuly**

Projet supervisé par **Elsa Dupraz**

## Remerciements

- **Elsa Dupraz** pour la supervision scientifique
- **Integrated DNA Technologies (IDT)** pour les données de modèle k4
- **DNArchive** pour les données de modèle k6
- La communauté de recherche en stockage ADN pour les avancées continues

---

## Ressources supplémentaires

### Documentation technique

- [BioPython Documentation](https://biopython.org/)
- [LDPC Codes Theory](https://en.wikipedia.org/wiki/Low-density_parity-check_code)
- [DNA Sequencing Technologies](https://www.genome.gov/about-genomics/fact-sheets/DNA-Sequencing-Fact-Sheet)

### Articles de référence

- Church, G. M., et al. (2012). "Next-Generation Digital Information Storage in DNA"
- Erlich, Y., & Zielinski, D. (2017). "DNA Fountain enables a robust and efficient storage architecture"
- Organick, L., et al. (2018). "Random access in large-scale DNA data storage"

### Projets connexes

- [Microsoft DNA Storage](https://www.microsoft.com/en-us/research/project/dna-storage/)
- [Catalog DNA](https://catalogdna.com/)
- [Twist Bioscience](https://www.twistbioscience.com/products/data-storage)

---

**DNA Storage Pipeline** - Transformer des mots en molécules, préserver l'information pour l'éternité.
