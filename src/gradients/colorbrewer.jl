#=

Apache-Style Software License for ColorBrewer software and ColorBrewer Color Schemes

Copyright (c) 2002 Cynthia Brewer, Mark Harrower, and The Pennsylvania State University.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

This text from my earlier Apache License Version 1.1 also remains in place for guidance on attribution and permissions:
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions as source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. The end-user documentation included with the redistribution, if any, must include the following acknowledgment:
"This product includes color specifications and designs developed by Cynthia Brewer (http://colorbrewer.org/)."
Alternately, this acknowledgment may appear in the software itself, if and wherever such third-party acknowledgments normally appear.
4. The name "ColorBrewer" must not be used to endorse or promote products derived from this software without prior written permission. For written permission, please contact Cynthia Brewer at cbrewer@psu.edu.
5. Products derived from this software may not be called "ColorBrewer", nor may "ColorBrewer" appear in their name, without prior written permission of Cynthia Brewer.

RGB values were taken from http://www.colorbrewer2.org
=#

register_color_library(:colorbrewer, ColorLibrary(Dict(:default => :sequential, :sequential => :Spectral, :divergent => :RdBu)))

register_gradient_colors(:Spectral, [
    RGB(158/255, 1/255, 66/255),
    RGB(213/255, 62/255, 79/255),
    RGB(244/255, 109/255, 67/255),
    RGB(253/255, 174/255, 97/255),
    RGB(254/255, 224/255, 139/255),
    RGB(255/255, 255/255, 191/255),
    RGB(230/255, 245/255, 152/255),
    RGB(171/255, 221/255, 164/255),
    RGB(102/255, 194/255, 165/255),
    RGB(50/255, 136/255, 189/255),
    RGB(94/255, 79/255, 162/255)], :colorbrewer)

register_gradient_colors(:RdYlGn, [
    RGB(165/255, 0/255, 38/255),
    RGB(215/255, 48/255, 39/255),
    RGB(244/255, 109/255, 67/255),
    RGB(253/255, 174/255, 97/255),
    RGB(254/255, 224/255, 139/255),
    RGB(255/255, 255/255, 191/255),
    RGB(217/255, 239/255, 139/255),
    RGB(166/255, 217/255, 106/255),
    RGB(102/255, 189/255, 99/255),
    RGB(26/255, 152/255, 80/255),
    RGB(0/255, 104/255, 55/255)], :colorbrewer)

register_gradient_colors(:RdBu, [
    RGB(103/255, 0/255, 31/255),
    RGB(178/255, 24/255, 43/255),
    RGB(214/255, 96/255, 77/255),
    RGB(244/255, 165/255, 130/255),
    RGB(253/255, 219/255, 199/255),
    RGB(247/255, 247/255, 247/255),
    RGB(209/255, 229/255, 240/255),
    RGB(146/255, 197/255, 222/255),
    RGB(67/255, 147/255, 195/255),
    RGB(33/255, 102/255, 172/255),
    RGB(5/255, 48/255, 97/255)], :colorbrewer)

register_gradient_colors(:PiYG, [
    RGB(142/255, 1/255, 82/255),
    RGB(197/255, 27/255, 125/255),
    RGB(222/255, 119/255, 174/255),
    RGB(241/255, 182/255, 218/255),
    RGB(253/255, 224/255, 239/255),
    RGB(247/255, 247/255, 247/255),
    RGB(230/255, 245/255, 208/255),
    RGB(184/255, 225/255, 134/255),
    RGB(127/255, 188/255, 65/255),
    RGB(77/255, 146/255, 33/255),
    RGB(39/255, 100/255, 25/255)], :colorbrewer)

register_gradient_colors(:PRGn, [
    RGB(64/255, 0/255, 75/255),
    RGB(118/255, 42/255, 131/255),
    RGB(153/255, 112/255, 171/255),
    RGB(194/255, 165/255, 207/255),
    RGB(231/255, 212/255, 232/255),
    RGB(247/255, 247/255, 247/255),
    RGB(217/255, 240/255, 211/255),
    RGB(166/255, 219/255, 160/255),
    RGB(90/255, 174/255, 97/255),
    RGB(27/255, 120/255, 55/255),
    RGB(0/255, 68/255, 27/255)], :colorbrewer)

register_gradient_colors(:RdYlBu, [
    RGB(165/255, 0/255, 38/255),
    RGB(215/255, 48/255, 39/255),
    RGB(244/255, 109/255, 67/255),
    RGB(253/255, 174/255, 97/255),
    RGB(254/255, 224/255, 144/255),
    RGB(255/255, 255/255, 191/255),
    RGB(224/255, 243/255, 248/255),
    RGB(171/255, 217/255, 233/255),
    RGB(116/255, 173/255, 209/255),
    RGB(69/255, 117/255, 180/255),
    RGB(49/255, 54/255, 149/255)], :colorbrewer)

register_gradient_colors(:BrBG, [
    RGB(84/255, 48/255, 5/255),
    RGB(140/255, 81/255, 10/255),
    RGB(191/255, 129/255, 45/255),
    RGB(223/255, 194/255, 125/255),
    RGB(246/255, 232/255, 195/255),
    RGB(245/255, 245/255, 245/255),
    RGB(199/255, 234/255, 229/255),
    RGB(128/255, 205/255, 193/255),
    RGB(53/255, 151/255, 143/255),
    RGB(1/255, 102/255, 94/255),
    RGB(0/255, 60/255, 48/255)], :colorbrewer)

register_gradient_colors(:RdGy, [
    RGB(103/255, 0/255, 31/255),
    RGB(178/255, 24/255, 43/255),
    RGB(214/255, 96/255, 77/255),
    RGB(244/255, 165/255, 130/255),
    RGB(253/255, 219/255, 199/255),
    RGB(255/255, 255/255, 255/255),
    RGB(224/255, 224/255, 224/255),
    RGB(186/255, 186/255, 186/255),
    RGB(135/255, 135/255, 135/255),
    RGB(77/255, 77/255, 77/255),
    RGB(26/255, 26/255, 26/255)], :colorbrewer)

register_gradient_colors(:PuOr, [
    RGB(127/255, 59/255, 8/255),
    RGB(179/255, 88/255, 6/255),
    RGB(224/255, 130/255, 20/255),
    RGB(253/255, 184/255, 99/255),
    RGB(254/255, 224/255, 182/255),
    RGB(247/255, 247/255, 247/255),
    RGB(216/255, 218/255, 235/255),
    RGB(178/255, 171/255, 210/255),
    RGB(128/255, 115/255, 172/255),
    RGB(84/255, 39/255, 136/255),
    RGB(45/255, 0/255, 75/255)], :colorbrewer)

# These are commented out as they constitute non-continuous gradients, and are thus used differently by Plots
# register_gradient_colors(:Set2, [
#     RGB(102/255, 194/255, 165/255),
#     RGB(252/255, 141/255, 98/255),
#     RGB(141/255, 160/255, 203/255),
#     RGB(231/255, 138/255, 195/255),
#     RGB(166/255, 216/255, 84/255),
#     RGB(255/255, 217/255, 47/255),
#     RGB(229/255, 196/255, 148/255),
#     RGB(179/255, 179/255, 179/255)], :colorbrewer)
#
# register_gradient_colors(:Accent, [
#     RGB(127/255, 201/255, 127/255),
#     RGB(190/255, 174/255, 212/255),
#     RGB(253/255, 192/255, 134/255),
#     RGB(255/255, 255/255, 153/255),
#     RGB(56/255, 108/255, 176/255),
#     RGB(240/255, 2/255, 127/255),
#     RGB(191/255, 91/255, 23/255),
#     RGB(102/255, 102/255, 102/255)], :colorbrewer)
#
# register_gradient_colors(:Set1, [
#     RGB(228/255, 26/255, 28/255),
#     RGB(55/255, 126/255, 184/255),
#     RGB(77/255, 175/255, 74/255),
#     RGB(152/255, 78/255, 163/255),
#     RGB(255/255, 127/255, 0/255),
#     RGB(255/255, 255/255, 51/255),
#     RGB(166/255, 86/255, 40/255),
#     RGB(247/255, 129/255, 191/255),
#     RGB(153/255, 153/255, 153/255)], :colorbrewer)
#
# register_gradient_colors(:Set3, [
#     RGB(141/255, 211/255, 199/255),
#     RGB(255/255, 255/255, 179/255),
#     RGB(190/255, 186/255, 218/255),
#     RGB(251/255, 128/255, 114/255),
#     RGB(128/255, 177/255, 211/255),
#     RGB(253/255, 180/255, 98/255),
#     RGB(179/255, 222/255, 105/255),
#     RGB(252/255, 205/255, 229/255),
#     RGB(217/255, 217/255, 217/255),
#     RGB(188/255, 128/255, 189/255),
#     RGB(204/255, 235/255, 197/255),
#     RGB(255/255, 237/255, 111/255)], :colorbrewer)
#
# register_gradient_colors(:Dark2, [
#     RGB(27/255, 158/255, 119/255),
#     RGB(217/255, 95/255, 2/255),
#     RGB(117/255, 112/255, 179/255),
#     RGB(231/255, 41/255, 138/255),
#     RGB(102/255, 166/255, 30/255),
#     RGB(230/255, 171/255, 2/255),
#     RGB(166/255, 118/255, 29/255),
#     RGB(102/255, 102/255, 102/255)], :colorbrewer)
#
# register_gradient_colors(:Paired, [
#     RGB(166/255, 206/255, 227/255),
#     RGB(31/255, 120/255, 180/255),
#     RGB(178/255, 223/255, 138/255),
#     RGB(51/255, 160/255, 44/255),
#     RGB(251/255, 154/255, 153/255),
#     RGB(227/255, 26/255, 28/255),
#     RGB(253/255, 191/255, 111/255),
#     RGB(255/255, 127/255, 0/255),
#     RGB(202/255, 178/255, 214/255),
#     RGB(106/255, 61/255, 154/255),
#     RGB(255/255, 255/255, 153/255),
#     RGB(177/255, 89/255, 40/255)], :colorbrewer)
#
# register_gradient_colors(:Pastel2, [
#     RGB(179/255, 226/255, 205/255),
#     RGB(253/255, 205/255, 172/255),
#     RGB(203/255, 213/255, 232/255),
#     RGB(244/255, 202/255, 228/255),
#     RGB(230/255, 245/255, 201/255),
#     RGB(255/255, 242/255, 174/255),
#     RGB(241/255, 226/255, 204/255),
#     RGB(204/255, 204/255, 204/255)], :colorbrewer)
#
# register_gradient_colors(:Pastel1, [
#     RGB(251/255, 180/255, 174/255),
#     RGB(179/255, 205/255, 227/255),
#     RGB(204/255, 235/255, 197/255),
#     RGB(222/255, 203/255, 228/255),
#     RGB(254/255, 217/255, 166/255),
#     RGB(255/255, 255/255, 204/255),
#     RGB(229/255, 216/255, 189/255),
#     RGB(253/255, 218/255, 236/255),
#     RGB(242/255, 242/255, 242/255)], :colorbrewer)

register_gradient_colors(:OrRd, [
    RGB(255/255, 247/255, 236/255),
    RGB(254/255, 232/255, 200/255),
    RGB(253/255, 212/255, 158/255),
    RGB(253/255, 187/255, 132/255),
    RGB(252/255, 141/255, 89/255),
    RGB(239/255, 101/255, 72/255),
    RGB(215/255, 48/255, 31/255),
    RGB(179/255, 0/255, 0/255),
    RGB(127/255, 0/255, 0/255)], :colorbrewer)

register_gradient_colors(:PuBu, [
    RGB(255/255, 247/255, 251/255),
    RGB(236/255, 231/255, 242/255),
    RGB(208/255, 209/255, 230/255),
    RGB(166/255, 189/255, 219/255),
    RGB(116/255, 169/255, 207/255),
    RGB(54/255, 144/255, 192/255),
    RGB(5/255, 112/255, 176/255),
    RGB(4/255, 90/255, 141/255),
    RGB(2/255, 56/255, 88/255)], :colorbrewer)

register_gradient_colors(:BuPu, [
    RGB(247/255, 252/255, 253/255),
    RGB(224/255, 236/255, 244/255),
    RGB(191/255, 211/255, 230/255),
    RGB(158/255, 188/255, 218/255),
    RGB(140/255, 150/255, 198/255),
    RGB(140/255, 107/255, 177/255),
    RGB(136/255, 65/255, 157/255),
    RGB(129/255, 15/255, 124/255),
    RGB(77/255, 0/255, 75/255)], :colorbrewer)

register_gradient_colors(:Oranges, [
    RGB(255/255, 245/255, 235/255),
    RGB(254/255, 230/255, 206/255),
    RGB(253/255, 208/255, 162/255),
    RGB(253/255, 174/255, 107/255),
    RGB(253/255, 141/255, 60/255),
    RGB(241/255, 105/255, 19/255),
    RGB(217/255, 72/255, 1/255),
    RGB(166/255, 54/255, 3/255),
    RGB(127/255, 39/255, 4/255)], :colorbrewer)

register_gradient_colors(:BuGn, [
    RGB(247/255, 252/255, 253/255),
    RGB(229/255, 245/255, 249/255),
    RGB(204/255, 236/255, 230/255),
    RGB(153/255, 216/255, 201/255),
    RGB(102/255, 194/255, 164/255),
    RGB(65/255, 174/255, 118/255),
    RGB(35/255, 139/255, 69/255),
    RGB(0/255, 109/255, 44/255),
    RGB(0/255, 68/255, 27/255)], :colorbrewer)

register_gradient_colors(:YlOrBr, [
    RGB(255/255, 255/255, 229/255),
    RGB(255/255, 247/255, 188/255),
    RGB(254/255, 227/255, 145/255),
    RGB(254/255, 196/255, 79/255),
    RGB(254/255, 153/255, 41/255),
    RGB(236/255, 112/255, 20/255),
    RGB(204/255, 76/255, 2/255),
    RGB(153/255, 52/255, 4/255),
    RGB(102/255, 37/255, 6/255)], :colorbrewer)

register_gradient_colors(:YlGn, [
    RGB(255/255, 255/255, 229/255),
    RGB(247/255, 252/255, 185/255),
    RGB(217/255, 240/255, 163/255),
    RGB(173/255, 221/255, 142/255),
    RGB(120/255, 198/255, 121/255),
    RGB(65/255, 171/255, 93/255),
    RGB(35/255, 132/255, 67/255),
    RGB(0/255, 104/255, 55/255),
    RGB(0/255, 69/255, 41/255)], :colorbrewer)

register_gradient_colors(:Reds, [
    RGB(255/255, 245/255, 240/255),
    RGB(254/255, 224/255, 210/255),
    RGB(252/255, 187/255, 161/255),
    RGB(252/255, 146/255, 114/255),
    RGB(251/255, 106/255, 74/255),
    RGB(239/255, 59/255, 44/255),
    RGB(203/255, 24/255, 29/255),
    RGB(165/255, 15/255, 21/255),
    RGB(103/255, 0/255, 13/255)], :colorbrewer)

register_gradient_colors(:RdPu, [
    RGB(255/255, 247/255, 243/255),
    RGB(253/255, 224/255, 221/255),
    RGB(252/255, 197/255, 192/255),
    RGB(250/255, 159/255, 181/255),
    RGB(247/255, 104/255, 161/255),
    RGB(221/255, 52/255, 151/255),
    RGB(174/255, 1/255, 126/255),
    RGB(122/255, 1/255, 119/255),
    RGB(73/255, 0/255, 106/255)], :colorbrewer)

register_gradient_colors(:Greens, [
    RGB(247/255, 252/255, 245/255),
    RGB(229/255, 245/255, 224/255),
    RGB(199/255, 233/255, 192/255),
    RGB(161/255, 217/255, 155/255),
    RGB(116/255, 196/255, 118/255),
    RGB(65/255, 171/255, 93/255),
    RGB(35/255, 139/255, 69/255),
    RGB(0/255, 109/255, 44/255),
    RGB(0/255, 68/255, 27/255)], :colorbrewer)

register_gradient_colors(:YlGnBu, [
    RGB(255/255, 255/255, 217/255),
    RGB(237/255, 248/255, 177/255),
    RGB(199/255, 233/255, 180/255),
    RGB(127/255, 205/255, 187/255),
    RGB(65/255, 182/255, 196/255),
    RGB(29/255, 145/255, 192/255),
    RGB(34/255, 94/255, 168/255),
    RGB(37/255, 52/255, 148/255),
    RGB(8/255, 29/255, 88/255)], :colorbrewer)

register_gradient_colors(:Purples, [
    RGB(252/255, 251/255, 253/255),
    RGB(239/255, 237/255, 245/255),
    RGB(218/255, 218/255, 235/255),
    RGB(188/255, 189/255, 220/255),
    RGB(158/255, 154/255, 200/255),
    RGB(128/255, 125/255, 186/255),
    RGB(106/255, 81/255, 163/255),
    RGB(84/255, 39/255, 143/255),
    RGB(63/255, 0/255, 125/255)], :colorbrewer)

register_gradient_colors(:GnBu, [
    RGB(247/255, 252/255, 240/255),
    RGB(224/255, 243/255, 219/255),
    RGB(204/255, 235/255, 197/255),
    RGB(168/255, 221/255, 181/255),
    RGB(123/255, 204/255, 196/255),
    RGB(78/255, 179/255, 211/255),
    RGB(43/255, 140/255, 190/255),
    RGB(8/255, 104/255, 172/255),
    RGB(8/255, 64/255, 129/255)], :colorbrewer)

register_gradient_colors(:Greys, [
    RGB(255/255, 255/255, 255/255),
    RGB(240/255, 240/255, 240/255),
    RGB(217/255, 217/255, 217/255),
    RGB(189/255, 189/255, 189/255),
    RGB(150/255, 150/255, 150/255),
    RGB(115/255, 115/255, 115/255),
    RGB(82/255, 82/255, 82/255),
    RGB(37/255, 37/255, 37/255),
    RGB(0/255, 0/255, 0/255)], :colorbrewer)

register_gradient_colors(:YlOrRd, [
    RGB(255/255, 255/255, 204/255),
    RGB(255/255, 237/255, 160/255),
    RGB(254/255, 217/255, 118/255),
    RGB(254/255, 178/255, 76/255),
    RGB(253/255, 141/255, 60/255),
    RGB(252/255, 78/255, 42/255),
    RGB(227/255, 26/255, 28/255),
    RGB(177/255, 0/255, 38/255)], :colorbrewer)

register_gradient_colors(:PuRd, [
    RGB(247/255, 244/255, 249/255),
    RGB(231/255, 225/255, 239/255),
    RGB(212/255, 185/255, 218/255),
    RGB(201/255, 148/255, 199/255),
    RGB(223/255, 101/255, 176/255),
    RGB(231/255, 41/255, 138/255),
    RGB(206/255, 18/255, 86/255),
    RGB(152/255, 0/255, 67/255),
    RGB(103/255, 0/255, 31/255)], :colorbrewer)

register_gradient_colors(:Blues, [
    RGB(247/255, 251/255, 255/255),
    RGB(222/255, 235/255, 247/255),
    RGB(198/255, 219/255, 239/255),
    RGB(158/255, 202/255, 225/255),
    RGB(107/255, 174/255, 214/255),
    RGB(66/255, 146/255, 198/255),
    RGB(33/255, 113/255, 181/255),
    RGB(8/255, 81/255, 156/255),
    RGB(8/255, 48/255, 107/255)], :colorbrewer)

register_gradient_colors(:PuBuGn, [
    RGB(255/255, 247/255, 251/255),
    RGB(236/255, 226/255, 240/255),
    RGB(208/255, 209/255, 230/255),
    RGB(166/255, 189/255, 219/255),
    RGB(103/255, 169/255, 207/255),
    RGB(54/255, 144/255, 192/255),
    RGB(2/255, 129/255, 138/255),
    RGB(1/255, 108/255, 89/255),
    RGB(1/255, 70/255, 54/255)], :colorbrewer)
