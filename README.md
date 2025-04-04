# Éléments finis: projet 2024-2025

# Groupe 139
- **Amory Geoffroy**
- **Desneux Tristan**

# Structure
Notre projet suit une structure plutôt classique: un dossier src avec le code source, un dossier data contenant divers maillage précalculés et un fichier defaultProblem.txt qui contient les paramètres physiques du problème d'élasticité linéaire. Un fichier CMakeLists.txt s'occupe de l'aspect compilation.

# Compilation
Pour compiler le projet, il est d'abord nécessaire d'ajouter un dossier glfw qui contient la librairie glfw et un dossier gmsh qui contient la librairie gmsh. À l'intérieur de ces dossiers, la structure est la même que pour les devoirs. Ensuite, une procédure classique avec un dossier build peut être utilisée pour compiler le programme.

# Exécution du programme
Après compilation, un fichier binaire myFem est crée, celui-ci permet d'exécuter le programme d'éléments finis. Plusieurs options (à ajouter à la commande dans le terminal) différents peuvent être utilisées pour diversifier l'exécution:
- l'option "-in *filepath*" permet de changer le maillage utilisé par le programme. Une grande variété de maillages sont disponibles et le programme s'adapte lui-même au types d'éléments utilisés (triangles ou quads)
- l'option "-p *filepath*" permet de changer les paramètres physiques du problème d'élasticité linéaire. Par contre le fichier fourni doit avoir strictement la même structure que le fichier defaultProblem.txt
- l'option "-gauss" permet de résoudre le problème avec un solveur classique qui effectue l'élimination gaussienne. Si cette option n'est pas sélectionnée, un solver band rcmk est utilisé.
- l'option "-axisym" résout le problème axysymmétrique.
- l'option "-dp" résout le problème en déformations planes. par défaut le problème est résolu en tensions planes pour mieux correspondre à notre application.

# Résultats
Lors de l'exécution du programme, un plot avec le champ de déformations s'affiche automatiquement, de plus les résultats numérique importants (min, max) sont données dans le terminal.

Il y a egalement des fichiers avec les contraines et les elongation sur chaque element qui se creent dans build. Pour obtenir la contraine maximale, executer givemaxcontr.py. 
