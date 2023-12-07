## 说明

### nanowires_fit_2_9_nm
nanowires_fit_2_9nm.m  是拟合主程序

程序 3-8 行是 Reference 数据 kpoints 和 energies.

如果只拟合最高价带的话，第6行代码为 "energies = energies(544:544,:)""， 以及 nanowires_valence.m 里的 N_E 变量改为 1.
如果拟合前10条最高价带，第6行代码为 "energies = energies(535:544,:)""， 以及 nanowires_valence.m 里的 N_E 变量改为 10.




### plot_nanowires_valence.m 
输入参数 LMN 可画出对应的能带。 生成的能带条数与 nanowires_valence.m 程序里面的 N_E 变量有关。