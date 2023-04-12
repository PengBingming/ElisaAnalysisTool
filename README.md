# ElisaAnalysisTool
A tool to analysis the Elisa data.

此应用已布置在 Hiplot平台：[ElisaAnalysis](https://hiplot.com.cn/cloud-tool/data-analysis/link/655)

功能：这是一个shiny-app，用于分析酶标仪的检测数据，app/app.R 为文件，file/ 文件夹中为参考数据。
安装相应包后即可使用。

使用：依据参考数据理解，其中前两列列名为 x 和 y 列，分别对应为 浓度 和 OD 值，其他列为所测样本 OD 值。

R包： xlsx readxl rstatix ggpubr ggplot2 shinydashboard shiny
注意：其中 xlsx 包在windows 可能无法使用，因此建议在 Linux 中使用。

