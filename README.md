# Overview

Seja bem vindo a este curso introdutório sobre análise de _variantes clínicas em DNA humano com dados NGS_. Durante este curso você irá aprender os conceitos básicos de análise e visualização de dados NGS e visualização em um ambiente LINUX/UNIX, os softwares mais usados e melhores práticas serão explanadas. Este curso é focado em alinhamento de dados para apenas DNA, chamada de variantes e visualização de dados.

Este workshop foi projetado para ser realizado em 4 encontros de 2 horas e meia e assume que todos tenham uma base sobre dados NGS e ambiente LINUX. Todos os materiais deste workshop são abertos e gratuitos, fique à vontade em usá-los. Todos os dados para os tutoriais estão disponíveis em nossa pasta do Dropbox [DropBox folder](https://www.dropbox.com/sh/4qkqch7gyt888h7/AABD_i9ShwryfAqGeJ0yqqF3a).

#Abstract

Avanços no Sequenciamento de Nova Geração (NGS) têm possibilitado oportunidades sem precedentes para a disciplina mineração de dados genéticos de pessoas a populações. A análise subsequente dos dados em busca de variantes genéticas que podem estar associadas com algum tipo de doença genética é um importante passo para novas alternativas de diagnóstico e tratamento personalizado de um paciente. Este workshop visa apresentar como realizar esta análise sobre os dados genéticos por meio da disciplina de bioinformática passando pelas etapas de alinhamento, chamada de variantes e anotação. O resultado final é um arquivo rico em informações que permitem aos médicos e pacientes conhecerem mais sobre seu DNA e auxiliá-los nas decisões relacionadas à sua saúde.


# Agenda

## Dia 1

Neste primeiro dia iremos aprender o básico com uma _introdução a Genética e a tecnologia NGS_ , isto inclui preparação de dados para a análise.  Você também irá aprender os primeiros passos com _shell LINUX/GNU_.  Vamos também aprender sobre o pipeline _Variant Calling_ e manipular os primeiros dados ligados a _Controle de Qualidade (QC)_ .

### Apresentação (15 min)
- Slides [[PDF](Course_Materials/presentation/SummerCourse_Apresentacao.pdf)]

### Introdução a Tecnologias NGS's (45 min)
Uma introdução rápida sobre as tecnologias NGS e alguns conceitos.
- [Slides](Course_Materials/intro-ngs/SummerCourse_NGSTechnologies.pdf)

### Introdução ao Linux/SHELL (50 min)
- Presentation [[PDF](Course_Materials/intro-linux/intro_Linux_mda14.pdf) | [OpenOffice](Course_Materials/intro-linux/intro_Linux_mda14.odp)]

### Controle de Qualidade para dados brutos de NGS (FASTQ) (30 min)
- [Slides](Course_Materials/quality_control/presentation/quality_control_presentation.pdf) - [Tutorial](Course_Materials/quality_control/tutorial/quality_control.html) - [Dados](https://www.dropbox.com/sh/4qkqch7gyt888h7/AAAqebBSC6JgDGq4emwNORCaa/quality_control)


## Dia 2

Neste segundo dia vamos focar em _Variant Calling_ e _Visualização dos dados_ com  alinhamento dos dados e análise das sequências.


### Alinhamento dos dados
- Slides [[PDF](Course_Materials/alignment/presentation/ngs-read-mapping-imedina-mda14.pdf) | - [Tutorial](Course_Materials/alignment/tutorial/example.html)

### Visualização dos dados NGS PT1

- [Slides](Course_Materials/visualization/presentation/2014-Cambridge_visualisation.pdf) || [Guia de Exemplo](Course_Materials/visualization/tutorial/000_example.html) || [Tutorial 1](Course_Materials/visualization/tutorial/010_example.html) || [Tutorial 2](Course_Materials/visualization/tutorial/020_example.html)

## Dia 3

O terceiro dia será voltado para a limpeza e recalibração dos dados de alinhamento e para a chamada de variantes gerando os arquivos VCF.

### Análise de variantes

- [Slides](Course_Materials/variant_calling/presentation/2014-Cambridge_variant_calling.pdf) || [Tutorial 1](Course_Materials/variant_calling/tutorial/010_example.html) || [Tutorial 2](Course_Materials/variant_calling/tutorial/020_example.html) || [Tutorial 3](Course_Materials/variant_calling/tutorial/030_example.html)
- 
### Análise de Cobertura

- [Slides](Course_Materials/association_studies/presentation/association_studies_presentation.pdf) - [Tutorial](Course_Materials/association_studies/tutorial/association_studies.html)

### Anotação de variantes

- [Slides](Course_Materials/variant_annotation/presentation/2014-Cambridge_variant_annotation.pdf) || [Annovar](Course_Materials/variant_annotation/tutorial/annovar.html) || [HPG-Variant](Course_Materials/variant_annotation/tutorial/hpg-variant.html)

## Dia 4

Neste último encontro o foco será na anotação das variantes com alguns bancos de dados e a detecção e filtragem de variantes relevantes ao caso clínico investigado.

### Priorização de variantes

- [Slides](Course_Materials/variant_prioritization/presentation/2014-Cambridge_variant_prioritization.pdf) || [Data](https://www.dropbox.com/sh/4qkqch7gyt888h7/AADPzrs9NGg0PjVqnwQocUJUa/annotation/hpg-variant/examples)

### Análises, Big Data. 

- [Slides](Course_Materials/variant_prioritization/presentation/2014-Cambridge_variant_prioritization.pdf)


### Para onde podemos seguir a partir de  agora ? 

- [Slides](Course_Materials/variant_prioritization/presentation/2014-Cambridge_variant_prioritization.pdf)


----

# Sobre

Este curso está sendo instruído pelos colaboradores da Genomika Diagnósticos. Você pode realizar qualquer pergunta para Marcel Caraciolo (marcel@genomika.com.br),  Rodrigo (rodrigo@genomika.com.br) e Filiphe Villar (filiphe@genomika.com.br).
