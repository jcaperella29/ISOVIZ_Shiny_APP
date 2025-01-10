# JCAP-Intron-Exon Visualization and Guide Table Using Isoviz

This app provides a user-friendly interface for the **Isoviz** R package, allowing users to easily create visualizations and tables for analyzing intron-exon structures and junction efficiency predictions. For detailed documentation on the Isoviz package, visit the [official GitHub repository](https://github.com/daklab/isoviz).

## Key Features

- **Junction to Isoform Map**: Generate and download a clear, visual map of exon-intron relationships in your gene of interest.
- **Guide Table Creation**: Produce a guide table with junction gRNA efficiency predictions (TIGER model) and export it as a CSV file.
- **Intuitive Interface**: Upload required files, customize settings, and generate outputs with ease.
- **Deployed Online**: Access the app directly at [ISOVIZ Shiny App Deployment](https://jcaperella1162.shinyapps.io/ISOVIZ_Shiny_APP/).

## How to Use the App

### Running the App
To use the app:
1. Access the deployed version at [ISOVIZ Shiny App Deployment](https://jcaperella1162.shinyapps.io/ISOVIZ_Shiny_APP/).
2. Alternatively, run it locally by placing the `app.R` file and `.css` file in the same folder and launching the `app.R` file in R or RStudio.

---

### Input Requirements

#### File Uploads
The following file inputs are required to generate outputs:
1. **Genome .psl File**: Upload your file containing expressed isoforms.
2. **Gene-Transcript Conversion .txt File**: Upload gene-to-transcript conversion data.
3. **Junction File (.junc.txt)**: Upload your LeafCutter junctions file.
4. **Intron Annotations .rda File**: Upload intron annotations.

#### Text Inputs
1. **Gene Name**: Enter the HGNC symbol of your gene of interest.
2. **Gene Ensembl ID**: Provide the Ensembl ID of your gene of interest.

#### Additional Inputs
1. **Junctions List (.txt)**: (Optional) Upload a list of junctions you want the guide table to focus on.
2. **Cell Type**: Specify the cell type of your samples (default: "custom").
3. **Minimum Junction Usage**: Set the minimum number of events required for a junction to be considered valid.

---

### Generating Outputs

#### Visualization
- **Generate Plot**: Click to create and display the Junction to Isoform map under the **Plot** tab. A notification will indicate when the process starts and completes. 
- **Download Plot**: Save the map as a `.png` file using the **Download Plot** button.

#### Guide Table
- **Generate Guide Table**: Click to create a guide table with gRNA efficiency predictions. This table appears under the **Guide Table** tab. Notifications will indicate the process's progress.
- **Download Table**: Save the table as a `.csv` file using the **Download Table** button.

---

### Notes
- When generating the plot or guide table, placeholder messages ("Your visualization will appear here" or "Your guide table will appear here") will be displayed until the content is generated. These messages will then appear below the results.

---

## Feedback and Contributions

We welcome feedback and suggestions to improve this app. Please feel free to open an issue or pull request on our [GitHub repository](https://github.com/your-repo-link). For questions or collaborations, contact us at jcaperella@gmail.com.
