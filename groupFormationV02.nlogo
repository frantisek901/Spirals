;; EXPERIMENTAL BRANCH!!!! EXPERIMENTAL FEATURES!!!!!!!!
;;
;; This version of our joint model goes further and study the differences in opinion distributions
;; We form two camps -- C1 and C2 -- they could have their mean opinion <-1, +1> and their sigma <0, 1>,
;; but we could also mirror them: C1-average = - C2-average; C1-sigma = C2-sigma,
;; we could also determine size of C1 <0, 1>, whize size of C2 = 1 - C1-size.
;; Here the average opinion holds for all opinions defining opinion position, in principle we might flip averages from opinion to opinion, but situation is still same -- one camp has different opinion from the another camp.
;; the homophily/heterophily is applied as follows:
;; Homophily:   we take the first 'N-agents' * 'C1-size' agents and for them we apply 'C1-average' and 'C1-sigma' when randomly drawing opinions,
;;              the result is that whole segment of small-world has opinions drawn according the same principle and parameters.
;; Heterophily: every agent has 'C1-size' chance to use 'C1-average' and 'C1-sigma' for random drawing of opinions,
;;              and 1 - 'C1-size' chance to use 'C2-average' and 'C2-sigma', the resul is well mixed public sphere.
;;
;; We apply mainly HK model, but Deffuant is still implemented in code and we look also at it in more than 1D and look how agents adapt in >1D opinion space and whether they form groups.
;;

;; Created:  2021-10-21 FranCesko
;; Edited:   2021-10-28 FranCesko
;; Encoding: windows-1250
;; NetLogo:  6.1.1
;;

;; IDEA: What about simply employ Spiral of Silence?
;;       Just simply -- general parameter on scale (0; 1> and probability of speaking her attitude/opinion,
;;       baseline is p==1, everybody speaks always, if p==0.5 so everybody has 0.5 probability to speak her opinion/attitude at given step,
;;       if succeeds - speaks in given step, if not - falls silent for the respective step.
;;       In HK mechanism, agent computes mean opinion of all speaking agents who are inside 'opinion boundary' (are not further than threshold).
;;       In Defuant, agent randomly takes one speaking agent inside the 'opinion boundary' and sets opinon as average of their opinions.
;; DONE!
;;
;; IDEA: Employ homophily/heterophily principle at model start.
;; In progress ...
;;
;; IDEA: Choose, how many opinions agents update: sometimes 1, 2, 3, 4 ...
;; DONE!
;;
;; IDEA: Compute clusters
;; DONE!
;;


;; WISHLIST:
;;     - differentiate between interpersonal communication and social media communication -- two overlapping networks with their own rules
;;     - how radicalization is possible? How polarization happens?
;;     - differential attraction
;;     - repulsion
;;     - media exposure will be crucial…we can ask abt opinion consistent content, opinion contrary, and “mainstream/mixed”…
;;       how to we conceptualize/model those in ABM? Is this too simplistic (eg, think of the different flavors of conservative media,
;;       ranging from CDU type media to extremist hate groups).
;;     - how to think about social media influencers (eg Trump before deplatforming)…
;;       is it possible to designate “superagents” who influence everyone sharing certain beliefs and see their effects…
;;       both reach everyone in a group and their opinions are very highly weighted (or people vary in how much they weight that opinion?
;;       Could estimate Twitter effect that way!  Perhaps one could even model how movement towards an opinion might influence the superagent
;;       to increase communication or change focus…
;;     - Give weights to opinions... Taken from media, or from interpersonal communication:
;;          -- agents pick opinion according the importance, and update importance according number of contacts regarding the opinion
;;
;;


;; Parameters:
;;   Small-world network (Watts-Strogatz)
;;   agents have more than 1 type of attitude
;;   opinion on scale <-1;+1>
;;   boundary -- defines range as fraction of maximum possible Eucleid distance in n-dimensional space, this maximum depends on number of opinions: sqrt(opinions * 4)
;;

;; TO-DO:
;; 1) stopping conditions -- DONE!
;;


extensions [nw]

turtles-own [Opinion-position Speak? Uncertainty Record Last-opinion]

globals [main-Record components positions]


;; Initialization and setup
to setup
  ca
  ask patches [set pcolor patch-color]

  if set-seed? [random-seed RS]

  nw:generate-watts-strogatz turtles links N-agents n-neis p-random [
    fd (max-pxcor - 1)
    set size (max-pxcor / 10)
    set Opinion-position get-opinion
    set Last-opinion Opinion-position
    set Record n-values record-length [0]
    set speak? speaking
    set Uncertainty get-uncertainty
    getColor
  ]

  ask patches [set pcolor patch-color]
  ask links [set hidden? TRUE]
  ask turtles [getPlace]
  set main-Record n-values record-length [0]

  reset-ticks
end

to-report get-opinion
  ;; original algorithm: n-values opinions [precision (1 - random-float 2) 3]
  let opinion n-values opinions [0]  ;; empty list -- we will use separate streams so we need define empty variable independently before
  ifelse homophily? [  ;; The first part is homophilic algorithm
    let breaking-point (C1-size * N-agents)
    let i 0
    while [i < opinions] [
      let o -2
      while [o < -1 or o > 1] [
        set o random-normal ifelse-value (who < breaking-point)[C1-average][C2-average] ifelse-value (who < breaking-point)[C1-sigma][C2-sigma]
      ]
      set opinion replace-item i opinion o
      set i (i + 1)
    ]
  ][  ;; The second part is heterophilic algorithm
    let C1? (C1-size > random-float 1)  ;; randomly ressolve whether agent is in C1 or C2
    let i 0
    while [i < opinions] [
      let o -2
      while [o < -1 or o > 1] [
        set o random-normal ifelse-value (C1?)[C1-average][ifelse-value (mirror-C1?) [0 - C1-average][C2-average]] ifelse-value (C1?)[C1-sigma][ifelse-value (mirror-C1?) [C1-sigma][C2-sigma]]
      ]
      set opinion replace-item i opinion o
      set i (i + 1)
    ]
  ]

  report opinion
end




to-report get-uncertainty
  let uValue 0
  if boundary-drawn = "constant" [set uValue boundary]
  if boundary-drawn = "uniform" [set uValue precision (random-float (2 * boundary)) 3]
  if boundary-drawn = "normal" [
    set uValue precision (random-normal boundary sigma) 3
    while [uValue < 0 or uValue > 1] [
      set uValue precision (random-normal boundary sigma) 3
    ]
  ]

  report uValue
end


to getPlace
  if X-opinion > opinions [set X-opinion 1]
  if Y-opinion > opinions [set Y-opinion 1]

  facexy ((item (X-opinion - 1) opinion-position) * max-pxcor) ((item (Y-opinion - 1) opinion-position) * max-pycor)
  set xcor (item (X-opinion - 1) opinion-position) * max-pxcor
  set ycor (item (Y-opinion - 1) opinion-position) * max-pycor
end

to getColor
  ifelse speak? [
    set color 15 + 4 * mean(opinion-position)
    set size (max-pxcor / 10)
  ][
    set color white
    set size 0
  ]
end

to-report speaking
  report p-speaking > random-float 1
end

to-report patch-color
  report 59.9 - (4.9 * (ln(1 + count turtles-here) / ln(N-agents)))
end

;; Main routine
to go

  if updating > opinions [set updating opinions]

  ask turtles [
    ;; speaking
    set speak? speaking
    set Last-opinion Opinion-position

    ;; Mechanism of opinion change
    ifelse model = "HK" [
      change-opinion-HK
    ][
      change-opinion-Deffuant
    ]

    ;; Getting place and coloring
    getColor
    getPlace
  ]

  ask patches [set pcolor patch-color]
  ask turtles [
    ;; we take 1 if opinion is same, we take 0 if opinion changes, then
    ;; we put on the start of the list Record, but we omit the last item from Record
    set Record fput ifelse-value (Last-opinion = Opinion-position) [1][0] but-last Record
  ]
  set main-Record fput precision (mean [mean Record] of turtles) 3 but-last main-Record

  tick

  if mean main-Record = 1 [last-steps stop]
end


to last-steps
  set components []
  set positions []
  while [count turtles > 0] [
    ask one-of turtles [
      let comp turtles with [opinion-position = [opinion-position] of myself]
      if count comp > smallest-component [
        set components fput count comp components
        set positions fput ([opinion-position] of one-of comp) positions
      ]
      ask comp [die]
    ]
  ]

  print components
  print positions
end


to change-opinion-HK

  let X-guys link-neighbors with [color != white]
  let lim-dist (Uncertainty * sqrt(opinions * 4))
  let influentials X-guys with [opinion-distance <= lim-dist]
  ;print lim-dist
  ;print sort [opinion-distance] of influentials

  if count influentials > 0 [
      let op-list shuffle (n-of updating (range opinions)) ;; here we draw a list of dimensions which we will update
      let steps updating
      let step 0
      while [step < steps] [
        let i item step op-list
        let their-op precision (mean [item i opinion-position] of influentials) 3
        set opinion-position replace-item i opinion-position their-op  ;; NOTE: H-K model really assumes that agent adopts immediatelly the 'consesual' position
        set step step + 1
      ]
    ;print step
  ]

 ; print opinion-position
 ; print [opinion-position] of influentials
 ; print [opinion-distance] of influentials
 ; print mean [item (X-opinion - 1) opinion-position] of influentials

end

to change-opinion-Deffuant
  let influentials link-neighbors with [color != white]
  set influentials influentials with [opinion-distance <= (Uncertainty * sqrt(opinions * 4))]

  if count influentials > 0 [

      let op-list shuffle n-of updating range opinions ;; here we draw a list of dimensions which we will update
      let partner one-of influentials
      let steps updating
      let step 0

      while [step < steps] [
        let myX item (item step op-list) opinion-position
        let herX item (item step op-list) [opinion-position] of partner
        let myU Uncertainty
        let herU [Uncertainty] of partner
        let Hij (min list (myX + myU) (herX + herU)) - (max list (myX - myU) (herX - herU))

        if Hij > herU [
          set opinion-position replace-item (item step op-list) opinion-position precision (myX + (mu * ((Hij / herU) - 1) * (herX - myX))) 3
          set Uncertainty precision (myU + (mu * ((Hij / herU) - 1) * (herU - myU))) 3
        ]
        set step step + 1
      ]
    ]
end


to-report opinion-distance
  let my opinion-position
  let her [opinion-position] of myself
  let steps length my
  let step 0
  let dist 0
  while [step < steps] [
    set dist dist + (item step my - item step her) ^ 2
    set step step + 1
  ]
  set dist sqrt dist
  ;print dist
  report dist
end
@#$#@#$#@
GRAPHICS-WINDOW
219
10
667
459
-1
-1
40.0
1
10
1
1
1
0
0
0
1
-5
5
-5
5
1
1
1
ticks
30.0

BUTTON
109
10
164
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
163
42
218
75
GO!
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
163
10
218
43
Step
go\n;ask turtles [\n;  set speak? TRUE\n;  getColor\n;]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
22
246
194
279
N-agents
N-agents
10
1000
400.0
10
1
NIL
HORIZONTAL

SLIDER
877
10
1049
43
n-neis
n-neis
1
150
50.0
1
1
NIL
HORIZONTAL

SLIDER
877
46
1049
79
p-random
p-random
0
0.5
0.02
0.01
1
NIL
HORIZONTAL

SLIDER
23
280
195
313
opinions
opinions
1
50
4.0
1
1
NIL
HORIZONTAL

BUTTON
27
469
82
502
getPlace
ask turtles [getPlace]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
4
10
96
55
model
model
"HK" "Deffuant"
0

SLIDER
23
312
195
345
p-speaking
p-speaking
0
1
1.0
0.01
1
NIL
HORIZONTAL

SLIDER
24
345
196
378
boundary
boundary
0.01
1
0.15
0.01
1
NIL
HORIZONTAL

INPUTBOX
676
10
748
70
RS
5.0
1
0
Number

SWITCH
749
10
859
43
set-seed?
set-seed?
0
1
-1000

BUTTON
84
469
139
502
getColor
ask turtles [\n  set speak? TRUE\n  getColor\n]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
140
471
195
504
inspect
inspect turtle 0\nask turtle 0 [print opinion-position]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
195
471
250
504
sizes
show sort remove-duplicates [count turtles-here] of turtles 
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
250
471
305
505
HIDE!
\nask turtles [set hidden? TRUE]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
303
473
358
506
SHOW!
\nask turtles [set hidden? FALSE]\nask turtles [set size (max-pxcor / 10)]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
4
88
96
133
boundary-drawn
boundary-drawn
"uniform" "constant"
0

SLIDER
24
410
196
443
sigma
sigma
0
1
0.05
0.01
1
NIL
HORIZONTAL

PLOT
1284
262
1484
382
Distribution of 'Uncertainty'
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.05 1 -16777216 true "" "histogram [Uncertainty] of turtles"

SLIDER
24
377
196
410
mu
mu
0.01
1
1.0
0.01
1
NIL
HORIZONTAL

BUTTON
358
473
449
506
avg. Uncertainty
show mean [Uncertainty] of turtles
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
445
473
507
506
Hide links
ask links [set hidden? TRUE]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
507
473
567
506
Show links
ask links [set hidden? FALSE]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
1055
10
1227
43
X-opinion
X-opinion
1
50
1.0
1
1
NIL
HORIZONTAL

SLIDER
1056
47
1228
80
Y-opinion
Y-opinion
1
50
2.0
1
1
NIL
HORIZONTAL

BUTTON
566
473
637
506
avg. Opinion
show mean [mean opinion-position] of turtles
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
674
258
1277
495
Developement of opinions
NIL
NIL
0.0
10.0
-0.01
0.01
true
true
"" ""
PENS
"Op01" 1.0 0 -16777216 true "" "plot mean [item 0 opinion-position] of turtles"
"Op02" 1.0 0 -7500403 true "" "if opinions >= 2 [plot mean [item 1 opinion-position] of turtles]"
"Op03" 1.0 0 -2674135 true "" "if opinions >= 3 [plot mean [item 2 opinion-position] of turtles]"
"Op04" 1.0 0 -955883 true "" "if opinions >= 4 [plot mean [item 3 opinion-position] of turtles]"
"Op05" 1.0 0 -6459832 true "" "if opinions >= 5 [plot mean [item 4 opinion-position] of turtles]"
"Op06" 1.0 0 -1184463 true "" "if opinions >= 6 [plot mean [item 5 opinion-position] of turtles]"
"Op07" 1.0 0 -10899396 true "" "if opinions >= 7 [plot mean [item 6 opinion-position] of turtles]"
"Op08" 1.0 0 -13840069 true "" "if opinions >= 8 [plot mean [item 7 opinion-position] of turtles]"
"Op09" 1.0 0 -14835848 true "" "if opinions >= 9 [plot mean [item 8 opinion-position] of turtles]"
"Op10" 1.0 0 -11221820 true "" "if opinions >= 10 [plot mean [item 9 opinion-position] of turtles]"
"Op11" 1.0 0 -13791810 true "" "if opinions >= 11 [plot mean [item 10 opinion-position] of turtles]"

SLIDER
4
55
96
88
updating
updating
1
50
2.0
1
1
NIL
HORIZONTAL

PLOT
674
107
1311
257
Stability of turtles (average)
NIL
NIL
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Turtles" 1.0 0 -16777216 true "" "plot mean [mean Record] of turtles"
"Main record" 1.0 0 -2674135 true "" "plot mean main-Record"

BUTTON
675
73
752
106
avg. Record
show mean [mean Record] of turtles
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
764
47
868
92
avg. Record
mean [mean Record] of turtles
6
1
11

SLIDER
1251
46
1423
79
record-length
record-length
10
100
100.0
1
1
NIL
HORIZONTAL

SLIDER
1252
13
1424
46
smallest-component
smallest-component
1
10
2.0
1
1
NIL
HORIZONTAL

SWITCH
96
75
218
108
homophily?
homophily?
1
1
-1000

SLIDER
96
140
218
173
C1-average
C1-average
-1
1
0.4
0.01
1
NIL
HORIZONTAL

SLIDER
96
173
218
206
C2-average
C2-average
-1
1
-0.3
0.01
1
NIL
HORIZONTAL

SLIDER
96
108
218
141
C1-size
C1-size
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
1
140
96
173
C1-sigma
C1-sigma
0
1
0.2
0.01
1
NIL
HORIZONTAL

SLIDER
1
173
96
206
C2-sigma
C2-sigma
0
1
0.15
0.01
1
NIL
HORIZONTAL

SWITCH
1338
86
1434
119
mirror-C1?
mirror-C1?
1
1
-1000

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)  

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="firstTry" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="5000"/>
    <metric>components</metric>
    <metric>positions</metric>
    <steppedValueSet variable="RS" first="1" step="1" last="10"/>
    <enumeratedValueSet variable="N-agents">
      <value value="129"/>
      <value value="257"/>
      <value value="513"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-random">
      <value value="0.05"/>
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-neis">
      <value value="4"/>
      <value value="16"/>
      <value value="64"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="opinions">
      <value value="1"/>
      <value value="2"/>
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="updating">
      <value value="1"/>
      <value value="8"/>
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary">
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-drawn">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;uniform&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-speaking">
      <value value="0.1"/>
      <value value="0.5"/>
      <value value="0.9"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="smallest-component">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X-opinion">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Y-opinion">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;HK&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-boundary">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sigma">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="set-seed?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mu">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="record-length">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
