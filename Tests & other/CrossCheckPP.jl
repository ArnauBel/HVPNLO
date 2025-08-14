
path_bdio_w = path_bdio_dict["local"]

# Initialize the LaTeX table string
latex_table = """
\\begin{table}[h]
\\centering
\\begin{tabular}{|c| c c c|}
\\hline
ens & a mPi & a mK & a fPi\\\\
\\hline
\\hline

"""

# Iterate over each ens.id
for (k,ens) in enumerate(ensInfo)
    row = ens.id * " & "

    mDict = BDIOread_mPP(path_bdio_w,ens.id)
    fDict = BDIOread_fPS(path_bdio_w,ens.id)

    difference = meson_ens[ens.id]["mPi"] - mDict["mPi"]; uwerr(difference)
    sigma = abs(value(difference))/err(difference)
    
    cell = "\\begin{tabular}{@{}c@{}} $(print_uwreal(meson_ens[ens.id]["mPi"])) \\\\ $(print_uwreal(mDict["mPi"])) \\\\ $(round(sigma,digits=2)) \\end{tabular}"
    row *= cell * " & "

    if ens.kappa_l == ens.kappa_s
        difference = meson_ens[ens.id]["mPi"] - mDict["mPi"]; uwerr(difference)
        sigma = abs(value(difference))/err(difference)
        
        cell = "\\begin{tabular}{@{}c@{}} $(print_uwreal(meson_ens[ens.id]["mPi"])) \\\\ $(print_uwreal(mDict["mPi"])) \\\\ $(round(sigma,digits=2)) \\end{tabular}"
        row *= cell * " & "
    else
        difference = meson_ens[ens.id]["mK"] - mDict["mK"]; uwerr(difference)
        sigma = abs(value(difference))/err(difference)
        
        cell = "\\begin{tabular}{@{}c@{}} $(print_uwreal(meson_ens[ens.id]["mK"])) \\\\ $(print_uwreal(mDict["mK"])) \\\\ $(round(sigma,digits=2)) \\end{tabular}"
        row *= cell * " & "
    end

    difference = meson_ens[ens.id]["fPi"] - fDict["fPi"]; uwerr(difference)
    sigma = abs(value(difference))/err(difference)

    cell = "\\begin{tabular}{@{}c@{}} $(print_uwreal(meson_ens[ens.id]["fPi"])) \\\\ $(print_uwreal(fDict["fPi"])) \\\\ $(round(sigma,digits=2)) \\end{tabular}"
    row *= cell

    row = row * " \\\\\n \\hline \n"
    latex_table *= row

    if k % 11 == 0 && k != length(ensInfo)
        latex_table *= """
        \\end{tabular}
        \\caption{Comparison with LD paper: mesons}
        \\end{table}
        """
        latex_table *= """
        \\begin{table}[h]
        \\centering
        \\begin{tabular}{|c| c c c|}
        \\hline
        ens & a mPi & a mK & a fPi\\\\
        \\hline
        \\hline

        """
    end
end

# Close the LaTeX table
latex_table *= """
\\end{tabular}
\\caption{Comparison with LD paper: mesons}
\\end{table}
"""


# Print or save the LaTeX table
println(latex_table)



