# Logboek Joris Deelstra

## 16-nov-2026

### Start-up

Vandaag is de start-up van deze module, en ook casus.

Het doel van deze casus is om een python module te maken, met erin classes heeft die tumorgroeimodellen kan maken en fitten.
Hierbij moet een uitgebreid markdown readme document met de volgende eisen:

1. Technische specificaties en afhankelijkheden
2. Bondige achtergrondinformatie over modellen en technieken
3. Referenties

Emiel heeft de github repo gemaakt, ik heb er de main.py en logboeken/ aan toegevoegd.

In de uitleg voor het project wordt er gerefereerd naar [DOI: 10.47611/jsrhs.v12i4.5202](https://doi.org/10.47611/jsrhs.v12i4.5202), echter is de DOI site op 16 nov niet bereikbaar. 
Ik heb artikel ook gevonden op de volgende site:
[ResearchGate.net](https://www.researchgate.net/publication/378286328_Treatment_Optimization_for_Tumor_Growth_by_Ordinary_Differential_Equations).

Het artikel gaat over hoe optimale behandelings tijdschemas kunnen worden gemaakt o.b.v. data van drie ODE modellen voor verschillende kanker soorten.
De conclussie is dat optimale tijdschemas per soort verschillend zijn, wat erop duid dat er geen universele is. Ook is er geen rekening gehouden met radiotherapie en bijwerkingen.
Wat voor ons interessant is is dat er veel modellen worden benoemd, en ook hoe modellen toegepast kunnen worden op kanker-groei; precies wat wij gaan doen. Ik neem aan dat we veel gaan refereren naar dit artikelen.

Nadat ik de casus documentatie heb doorgelezen heb ik samen met Emiel gezeten voor een taakverdeling.
Wat daaruit kwam, is dat we het afwisselend per workshop de bijbehorende modellen zullen toevoegen aan één python class.

Ik heb nu de basis class gemaakt. Gezien dat binnen twee minuutjes klaar is heb ik een plot-functie gemaakt die de data waarmee de class wordt ingeladen plot.
Emiel heeft dummy data gegenereerd.

### Data-Analyse

De data structuur is als volgt:
- Subject_ID (str) Group (str) Time_Days (int) Tumor_Volume_mm3 (float) Body_Weight_g (Float)

* Er zijn 30 Subjects
* Er zijn 3 groepen (control, treatment low, treatment high) met elk 10 subjects
* De tijdspan is 28 dagen, elke 3 dagen wordt er een sample genomen

Wat uit de data-analyse komt, waarin ik tumor volume tegen tijd (28 dagen) heb geplot, is als volgt:
1. Control neemt van ~100 naar ~1000-1400 toe, met 1 uitzondering waar op dag 15 een sample op 1400 zit, en een Subject waarbij de tumor langzaam afneemt van grootte.
2. Treatment Low neemt van ~100 naar ~500-700 toe, met 1 uitzondering op dag 17 waar een sample ~0 is, en een subject waarbij de tumor stopt met groeien op ~300.
3. Treatment High neemt ~niet toe, met een uitzondering op dag 10 waar een sample ~900 volume wordt, en van een subject die van ~100 naar ~1800 volume gaat.

Alhoewel het niet zo handig is voor sommige modellen dat er outliers zijn, gezien we de modellen moeten fitten, lijkt het ons wel interessant om de data zo te houden om te kijken hoeveel invloed het heeft op de modellen.
Na overleg met Emiel willen we misschien wel een tweede dataset aanmaken maar dan zonder de outliers.

### Workshop 2
De tweede workshop ging over ODE's.
Met ODE's kan je als je een twee dimensionale dataset hebt en per x stap meerdere y punten hebt het gemiddelde berekenen per x stap en dan een trendline kan maken in de grafiek.
Na de workshop ben ik bezig gegaan met het Euler model te implementeren in de class.
Het Euler model berekent telkens de afgeleide aan het begin van een x-stap (b.v. van x=0 naar x=1), en maakt op basis daarvan een y-stap.

Voor onze tumor data set werkt dit prima. De lijn komt relatief overeen met de data.
Ik heb het nog niet geprobeerd zonder de outliers eruit te filteren / n nieuwe dataset te genereren zonder outliers.

Bij de corona data-set gaf numpy divide by zero errors, en was de reference 'Recovered' niet geplot.
Dit komt omdat de regions als startpunt allemaal 0 hebben; hierdoor wordt r oneindig.
Deze error heb ik 'opgelost' door een if statement toe te voegen waarbij als het 0 is het wordt vervangen door 0.1.

Het Euler model-method werkt overigens ook niet bij deze dataset. Dit komt omdat de statistische model linear toe of af neemt, en niet trendlines kan maken die een andere vorm dan exponentieel is. Dit is vooral te zien bij Infected per Region. ![image](images/Euler_plot_sirdata_infected_per_region.png)

## 18-nov-2026

Vandaag ben ik bezig geweest met de andere ODE methodes toe te voegen.
Gezien ze op dezelfde principes werken als Euler heb ik een prepare functie aangemaakt, en dan per methode een aparte functie heb gemaakt die afhankelijk is van de prepare functie.

## 19-nov-2026

Ik ging vandaag eerst bezig met opzetten van de README. Ik heb het skelet gemaakt en ben begonnen met de hoofdstukken 1, 3 en referenties.
Daarna ging ik bezig met de wiskundige modellen te implementeren, die in brightspace staan.
Omdat ik en Emiel graag alle 11 wiskundige modellen willen implementeren in de jupyter notebook en gezien we niet genoeg tijd hebben om ze handmatig te implementeren, hebben we besloten om de wiskundige formules te laten schrijven door AI.

---

### **AI disclosure statement**

**Vraag aan Claude opus 5.5 (AI):**
Ik heb gevraagd om de groeimodellen uit de TODO in `mathematical_methods` toe te voegen. De functie `prepare` in `ODE_methods` moest die modellen gaan gebruiken in plaats van de vaste exponentiële groei. Daarna heb ik gevraagd om alle comments naar het Engels te vertalen.

**Wat er is gedaan:**
- `mathematical_methods` geeft nu voor elk van de 11 modellen de differentiaalvergelijking `f(t, V)` terug. De modellen zijn lineair, exponentieel, Mendelsohn, exponentieel afvlakkend, logistisch, Montroll, Allee, lineair gelimiteerd, oppervlakte-gelimiteerd, Von Bertalanffy en Gompertz.
- De nieuwe functie `estimate_params` schat startparameters uit de data. Vmax en Vmin volgen uit de grootste en kleinste meetwaarde. De groeisnelheid c wordt zo gekozen dat het model halverwege even snel stijgt als de data gemiddeld doet.
- `ODE_methods` heeft twee nieuwe argumenten gekregen. Met `growth_model` kies je het model en met `params` geef je zelf parameters op. `prepare` haalt de vergelijking nu op via `mathematical_methods`.
- Het script is getest door AI met alle modellen in combinatie met Euler, Heun en Runge-Kutta. De comments zijn daarna naar het Engels vertaald.

**Waarom zo:**
- Door de vergelijking uit `prepare` te halen, werkt elk model met elke solver. We hoeven dan niets dubbel te schrijven.
- Het exponentiële model is de standaard en gebruikt nog dezelfde schatting als eerst. Bestaande code geeft daardoor dezelfde uitkomst.
- De parameterschatting is bewust eenvoudig gehouden. Een goede fit hoort bij Workshop 3 (optimalisatie). De schatting is daar een startpunt voor, want de meeste modellen passen nog niet goed op de data. Emiel gaat hiet verder mee bezig.
- De modelnamen zijn Nederlands gebleven, zodat ze overeenkomen met de namen uit de opdracht.

---

## 22-nov-2026

### Workshop 3 & 4: optimalisatie en modelvergelijking

Vandaag heb ik `plot_optimization` aangepast zodat alle wiskundige modellen uit `mathematical_methods` erop gefit kunnen worden, in plaats van alleen het lineaire model.
Dit was nodig omdat `model_comparison` alle modellen achter elkaar moet fitten en vergelijken.

Wat ik heb aangepast:
1. Ik heb een hulpfunctie `_solve_ode` toegevoegd die de differentiaalvergelijking oplost met Euler, Heun of Runge-Kutta.
2. `plot_optimization` neemt nu een `model` (bijv. 'gompertz') en een `method` (de ODE-solver) als input. Per behandelgroep wordt het model opgelost en de MSE berekend tegen alle meetpunten.
3. De startwaarden komen uit `estimate_params`, berekend op het gemiddelde volume per tijdstip. Het startvolume is het gemiddelde op de eerste dag.
4. De random search uit Workshop 3 is gebleven, maar de stapgrootte is nu relatief per parameter. Met een vaste sigma van 0.01 werkte het niet, omdat `c` rond de 0.05 ligt en `Vmax` rond de 1000. Er is ook een maximum aantal pogingen toegevoegd zodat hij niet eindeloos doorloopt.
5. Er wordt alleen nog een figuur gemaakt als `plot=True`, anders kregen we bij de modelvergelijking 11 lege figuren.
6. De functie geeft nu een dict terug met de parameters en MSE per groep, de totale SSE, het aantal datapunten (n) en het aantal parameters (k).

Ook heb ik dus de functie `model_comparison` gemaakt. 
De functie loopt over elke type wiskundige model, runt optimization voor dat model, berekent de AIC AICc en BIC voor dat model en geeft dan een df waarin alle modellen worden vergelijkt.

Echter geeft het momenteel wel errors:

```error
C:\Users\Joris\OneDrive - Hanze\Bio-informatica\J3\P2\Casus A\case1_modelleren\main.py:186: RuntimeWarning: overflow encountered in scalar power
  "mendelsohn": lambda t, V: c * pos(V) ** d,
C:\Users\Joris\OneDrive - Hanze\Bio-informatica\J3\P2\Casus A\case1_modelleren\main.py:416: RuntimeWarning: overflow encountered in square
  mse = np.mean((y - y_model) ** 2)
C:\Users\Joris\OneDrive - Hanze\Bio-informatica\J3\P2\Casus A\case1_modelleren\main.py:190: RuntimeWarning: overflow encountered in scalar multiply
  "allee": lambda t, V: c * (V - Vmin) * (Vmax - V),
C:\Users\Joris\OneDrive - Hanze\Bio-informatica\J3\P2\Casus A\case1_modelleren\main.py:416: RuntimeWarning: overflow encountered in square
  mse = np.mean((y - y_model) ** 2)
C:\Users\Joris\OneDrive - Hanze\Bio-informatica\J3\P2\Casus A\case1_modelleren\main.py:376: RuntimeWarning: overflow encountered in scalar multiply
  ys[i + 1] = V + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
```

Gezien het nu al laat is, ga ik er morgen mee bezig

## 23-nov-2026

Vandaag ga ik bezig met de puntjes op de i te zetten voor het project, beginnend met dat de code momenteel overflow errors geeft.

### Overflow errors:

De overflows komen waarschijnlijk uit `optimization`, waar random.gauss instaat wat ruis toevoegd aan de parameters. Daardoor worden af en toe combinations geprobeerd waarbij de oplossing van de ODE naar oneindig gaat.
Ik ga paar vangnetten toevoegen in de code, om dit soort van errors te voorkomen.

Verassend genoeg was een kleine fix genoeg voor de meerdere errors.
Ik heb in `_solve_ode` een check toegevoegd wat kijkt of vars np.infinite worden of groter dan een belachelijk groot getal.
Deze fix heeft in een keer de 5 verschillende errors opgelost, chill.

### Installation

Ik heb een requirements.txt toegevoegd, en een installatie guide toegevoegd aan de README.

### Jupyter notebook bijwerken:

