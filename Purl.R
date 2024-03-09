doc = 2

library(knitr)
purl(input = "2A-BayesianLogic.qmd", 
     output = "classroom/1.1-Introduction/Section2InR.R",        
            documentation = doc)

purl(input = "2B-PosteriorEstimation.qmd", 
     output = "classroom/1.2-PosteriorEstimation/Section3InR.R",        
     documentation = doc)

purl(input = "2C-LMM.qmd", 
     output = "classroom/1.3-LMM/Section4InR.R",        
     documentation = doc)

purl(input = "3A-ModelSelection.qmd", 
     output = "classroom/2.1-BayesianModelSelection/Section5InR.R",        
     documentation = doc)

purl(input = "3B-Workflow.qmd", 
     output = "classroom/2.2-Workflow/Section6InR.R",        
     documentation = doc)

purl(input = "3C-GLMM.qmd", 
     output = "classroom/2.3-GLMM/Section7InR.R",        
     documentation = doc)

purl(input = "4A-ErrorInVariable.qmd", 
     output = "classroom/3.1-ErrorInVariableModels/Section8InR.R",        
     documentation = doc)

purl(input = "4B-Occupancy.qmd", 
     output = "classroom/3.2-OccupancyModels/Section9InR.R",        
     documentation = doc)

purl(input = "4C-StateSpaceModels.qmd", 
     output = "classroom/3.3-State-space/Section10InR.R",        
     documentation = doc)

purl(input = "4D-AutoregressiveModels.qmd", 
     output = "classroom/3.4-Autoregressive/Section11InR.R",        
     documentation = doc)

purl(input = "4E-IntegratedModels.qmd", 
     output = "classroom/3.5-IntegratedModels/Section12InR.R",        
     documentation = doc)

purl(input = "4F-BayesianSEMs.qmd", 
     output = "classroom/3.6-CausalModels/Section13InR.R",        
     documentation = doc)

purl(input = "4G-ProcessBased.qmd", 
     output = "classroom/3.7-Process-based/Section14InR.R",        
     documentation = doc)

purl(input = "4H-ApproximateBayesian.qmd", 
     output = "classroom/3.8-ApproximateBayesian/Section15InR.R",        
     documentation = doc)




#files = list.files(pattern = "*.qmd")
#for(i in 1:length(files)) knitr::purl(input = files[[i]], 
#            documentation = doc)
