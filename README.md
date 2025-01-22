# JCAP - Intron-Exon Visualization and Guide Table Using Isoviz

This app provides a user-friendly interface for the **Isoviz** R package, enabling users to create visualizations and tables for analyzing intron-exon structures and predicting junction efficiency. For detailed documentation on the Isoviz package, visit the [official GitHub repository](https://github.com/daklab/isoviz).

---

## 🚀 Key Features

- **Junction to Isoform Map**: Create and download clear, interactive visualizations of exon-intron relationships for your gene of interest.
- **Guide Table Creation**: Generate a guide table with junction gRNA efficiency predictions (using the TIGER model) and export it as a CSV file.
- **User-Friendly Interface**: Easily upload required files, customize analysis parameters, and generate outputs.

---

## 🛠 Current Status: Local Deployment Only

⚠️ **Note**: The app is currently available for local deployment only. Due to resource-intensive processes, it cannot be hosted on the free tier of shinyapps.io. Users are encouraged to run the app locally for now.

---

## 🛠 How to Use the App

### Running the App Locally
1. Clone the repository or download the app files.
2. Ensure the `app.R` file and all required files (e.g., `.css` and `.rda`) are in the same directory.
3. Run the app locally in R or RStudio:
   ```R
   library(shiny)
   runApp("path/to/app.R")

👾💻🥷 Got it! I'll update the README to indicate that the app currently works for local deployment only. Here's the updated version:

📝 Updated README
markdown
Copy
Edit
# JCAP - Intron-Exon Visualization and Guide Table Using Isoviz

This app provides a user-friendly interface for the **Isoviz** R package, enabling users to create visualizations and tables for analyzing intron-exon structures and predicting junction efficiency. For detailed documentation on the Isoviz package, visit the [official GitHub repository](https://github.com/daklab/isoviz).

---

## 🚀 Key Features

- **Junction to Isoform Map**: Create and download clear, interactive visualizations of exon-intron relationships for your gene of interest.
- **Guide Table Creation**: Generate a guide table with junction gRNA efficiency predictions (using the TIGER model) and export it as a CSV file.
- **User-Friendly Interface**: Easily upload required files, customize analysis parameters, and generate outputs.

---

## 🛠 Current Status: Local Deployment Only

⚠️ **Note**: The app is currently available for local deployment only. Due to resource-intensive processes, it cannot be hosted on the free tier of shinyapps.io. Users are encouraged to run the app locally for now.

---

## 🛠 How to Use the App

### Running the App Locally
1. Clone the repository or download the app files.
2. Ensure the `app.R` file and all required files (e.g., `.css` and `.rda`) are in the same directory.
3. Run the app locally in R or RStudio:
   ```R
   library(shiny)
   runApp("path/to/app.R")
📂 Input Requirements
Required File Uploads
Genome .psl File: A file containing expressed isoforms for the analysis.
Gene-Transcript Conversion .txt File: A file mapping genes to transcripts.
Junction File (.junc.txt): A LeafCutter junctions file.
Intron Annotations .rda File: An R data file containing intron annotations.
Additional Input Fields
Gene Name: The HGNC symbol of the gene to analyze (e.g., "RBFOX2").
Gene Ensembl ID: The Ensembl ID of the gene to analyze (e.g., "ENSG00000100320").
Junctions List (.txt): (Optional) A text file listing junctions to focus on.
Cell Type: Specify the cell type of your data (default: "Custom").
Minimum Junction Usage: Set the minimum number of events required for a junction to be considered (default: 5).
📊 Generating Outputs
Visualization
Generate Plot: Click the "Generate Plot" button to create the Junction to Isoform map, which will appear under the Plot tab.
Download Plot: Save the visualization as a .png file using the "Download Plot" button.
Guide Table
Generate Guide Table: Click the "Generate Guide Table" button to create a guide table with gRNA efficiency predictions. The table will appear under the Guide Table tab.
Download Table: Save the guide table as a .csv file using the "Download Table" button.
⚠️ Notes
Placeholder messages ("Your visualization will appear here" or "Your guide table will appear here") will display until the output is generated.
Long-running processes may take a few minutes. Please be patient.
🤝 Feedback and Contributions
We welcome feedback, suggestions, and contributions!

To report issues or suggest features, open an issue or pull request on our GitHub repository.
For questions or collaboration inquiries, contact us at jcaperella@gmail.com.
