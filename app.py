import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import re

@st.cache_data
def load_data(file_path, index_col = 0):
    return pd.read_csv(file_path, index_col = index_col)


def get_pval_stars(p):
    if pd.isna(p):
        return "n.s."
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "n.s."

uploaded_file = st.file_uploader("Upload an Excel file", type=["xlsx", "xls"])


if uploaded_file:
    #read excel file in and produce protein search 
    df = pd.read_excel(uploaded_file)
    st.write("Preview of uploaded data:", df.head())
    st.sidebar.header("Search For A Protein")

    #extract abundance column names and p-value column names
    abundance_cols = [col for col in df.columns if col.lower().startswith("abundances (grouped):")]
    pvalue_cols = [col for col in df.columns if col.lower().startswith("abundance ratio p-value:")]

    st.write("Detected Abundance Columns:", abundance_cols)
    st.write("Detected P-value Columns:", pvalue_cols)
    if not abundance_cols:
        st.error("No 'grouped abundance' columns found in the dataset. Please check your file.")
    else:
        protein_col = 'Accession'
        search_term = st.sidebar.text_input("Search:")

        if search_term:
            protein_data = df[df[protein_col] == search_term][abundance_cols].T
            protein_data.index = [col.replace("Abundances (Grouped): ", "") for col in protein_data.index]
            conditions = list(set([col.split(": ")[1] for col in abundance_cols]))

            if protein_data.empty:
                st.error(f"No data found for protein: {search_term}")
            else:
                protein_data.columns = ["Abundance"]
                protein_data.index.name = "Condition"
                protein_data.reset_index(inplace=True)

                pvalue_map = {}
                for col in pvalue_cols:
                    match = re.search(r"Abundance Ratio P-Value: \(([^)]+)\) / \(([^)]+)\)", col)
                    if match:
                        cond_A, cond_B = match.groups()
                        pvalue_map[(cond_A, cond_B)] = df[df[protein_col] == search_term][col].values[0]

                #Plot
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.bar(protein_data["Condition"], protein_data["Abundance"])
                ax.set_xlabel("Group")
                ax.set_ylabel("Grouped Abundances")
                ax.set_title(search_term)
                plt.xticks(rotation=45, ha="right")

                max_y = protein_data["Abundance"].max()
                offset = max_y * 0.05  


                num_comparisons = 0 
                #iterate through dictionary and draw brackets
                for (cond_A, cond_B), p_val in pvalue_map.items():
                    if cond_A in protein_data["Condition"].values and cond_B in protein_data["Condition"].values:
                        star_label = get_pval_stars(p_val)
                        if star_label == "n.s.":
                            continue
                        
                        x1 = protein_data[protein_data["Condition"] == cond_A].index[0]
                        x2 = protein_data[protein_data["Condition"] == cond_B].index[0]

                        # adjust bracket by offset to prevent overlap
                        y = max_y + offset + (num_comparisons * 0.1 * max_y)

                        num_comparisons += 1

                        # Draw bracket
                        ax.plot([x1, x1, x2, x2], [y, y + 0.02, y + 0.02, y], lw=1.5, color="black")

                        # add stars based on significance level
                        
                        ax.text((x1 + x2) / 2, y + 0.02, star_label, ha='center', va='bottom', fontsize=12, color="black")

    

                st.pyplot(fig)
        
        




