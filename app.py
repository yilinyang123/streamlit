import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
adata = sc.read("thymic_APCs_processed_00.h5ad")
gene_list = list(adata.var_names)
st.set_option('deprecation.showPyplotGlobalUse', False)
st.title('thymic_APCs_processed_00')

#user_input = st.text_input("gene", "Cd83")
user_input = st.multiselect("gene",gene_list,['Cd83'])
user_input = user_input[0]
coord = adata.obsm['X_umap']
cluster = pd.Categorical(adata.obs['leiden']).astype(int)
a = adata[: , user_input].X.toarray().astype('float')

@st.cache(allow_output_mutation=True)
def get_fig():
	fig,(ax1, ax2) = plt.subplots(1,2,figsize=(15,6))
	ax1.scatter(coord[:, 0], coord[:, 1], c = cluster, cmap = 'tab20', s = 2,alpha=0.65)
	sc = ax2.scatter(coord[:, 0], coord[:, 1], c = a, cmap = 'viridis', s = 2,alpha=0.65)
	fig.colorbar(sc, ax=ax2)
	ax1.axis('off')
	ax2.axis('off')
	ax1.set_title("Cluster")
	ax2.set_title(user_input)
	return fig

st.pyplot(get_fig())
#st.write(get_fig())
#st.altair_chart(fig)
