# EEG Classification for Motor Imagery Tasks  

This repository contains MATLAB code developed as part of my Master's thesis:  
**"Exploring Motor Imagery-Induced Neural Activation and Corticospinal Excitability in Healthy Adults Using EEG and TMS"**  
(Memorial University of Newfoundland, 2025).  

📄 Full thesis report: [Thesis_MonaHejazi.pdf](./Thesis_MonaHejazi.pdf)  

---

## 📌 Project Overview  
This project investigates whether EEG can be used to classify **motor imagery (MI)** tasks targeting **fine motor control** and **motor coordination**.  

Two feature extraction strategies were explored:  
- **Filter Bank Common Spatial Patterns (FBCSP)**  
- **Functional Connectivity (FC)** (correlation, coherence, PLV)  

Multiple classifiers were tested (LDA, SVM, Decision Trees) to evaluate performance across different pipelines.  

---

## ⚙️ Repository Structure  

```
EEG_MI_Classification/
│── README.md                     # Project documentation  
│── Thesis_MonaHejazi.pdf          # Full thesis report  
│── .gitignore                     # Git ignore file  
│
├── Preprocessing/
│   └── PreProcessing.m            # EEG preprocessing (downsampling, notch filter, bandpass filtering)  
│
├── FBCSP/                         # Feature extraction & classification using CSP/FBCSP  
│   ├── CSP_classification_5b_22ch.m  
│   ├── CSP_classification_5b_64ch.m  
│   ├── CSP_classification_9b_22ch.m  
│   └── CSP_classification_9b_64ch.m  
│
├── FunctionalConnectivity/        # Functional connectivity features & classification  
│   ├── Features/
│   │   ├── FC_allfeatures_5b_22ch.m  
│   │   ├── FC_allfeatures_5b_64ch.m  
│   │   ├── FC_allfeatures_9b_22ch.m  
│   │   └── FC_allfeatures_9b_64ch.m  
│   │
│   └── Classification/
│       ├── FC_classification_5b_22ch.m  
│       ├── FC_classification_5b_64ch.m  
│       ├── FC_classification_9b_22ch.m  
│       └── FC_classification_9b_64ch.m  
│
└── Statistical Analysis/
    └── AllPipelines_StatisticalAnalysis.m   # Statistical comparison of pipelines and visualization  
```

---

## 📊 Key Findings  
- **Preprocessing:** Downsampling, notch filter (60 Hz), and bandpass filtering were applied.  
- **FBCSP:** CSP-based features extracted across 5-band and 9-band configurations with 22 vs. 64 channels.  
- **Functional Connectivity:** Features derived from correlation, coherence, and PLV.  
- **Classification:** SVM, LDA, and Decision Trees evaluated.  
- **Best pipeline:** PLV-based functional connectivity achieved the highest accuracies (up to **80–83%** for some participants).  

---

## 📑 Citation  
If you use this code or report, please cite:  

Hejazi, Mona. *Exploring Motor Imagery-Induced Neural Activation and Corticospinal Excitability in Healthy Adults Using EEG and TMS.*  
Master’s Thesis, Memorial University of Newfoundland, July 2025.  

---

## 📬 Contact  
For questions or collaborations, feel free to reach out:  
📧 mhejazi@mun.ca 
🔗 [LinkedIn: https://www.linkedin.com/in/mona-hejazi-685592a5/]  
