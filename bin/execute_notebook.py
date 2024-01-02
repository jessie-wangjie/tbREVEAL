import papermill as pm

def create_report():
    pm.execute_notebook(
    '/data/tbHCA/bin/report_generation.ipynb',
    '/data/tbHCA/bin/papermill_test.ipynb',
    parameters = {'results_file':'/data/CM_HC_TB000200e_20231017_analysis/work/e9/db124346eb78d2ec28e2ee68e528b5/CM_HC_TB000200e_20231017_results.xlsx'}
)
    
if __name__ == "__main__":
    create_report()