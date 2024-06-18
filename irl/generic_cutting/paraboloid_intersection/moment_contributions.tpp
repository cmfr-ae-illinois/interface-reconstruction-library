// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2022 Fabien Evrard <fa.evrard@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_
#define IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_

#include <float.h>
#include <cassert>
#include <cmath>

#include "external/NumericalIntegration/NumericalIntegration.h"
#include "irl/data_structures/small_vector.h"
#include "irl/data_structures/stack_vector.h"
#include "irl/generic_cutting/half_edge_cutting/half_edge_cutting_helpers.h"
#include "irl/geometry/general/normal.h"
#include "irl/geometry/general/pt.h"
#include "irl/geometry/general/reference_frame.h"
#include "irl/geometry/general/rotations.h"
#include "irl/geometry/general/scalar_with_gradient.h"
#include "irl/geometry/general/unit_quaternion.h"
#include "irl/helpers/mymath.h"
#include "irl/paraboloid_reconstruction/paraboloid.h"

namespace IRL {

template <class ScalarType, UnsignedIndex_t Order>
const std::array<ScalarType, Order> AbscissaeGauss(void) {
  if constexpr (Order == 5) {
    return std::array<ScalarType, Order>(
        {ScalarType(-0.9061798459386639927976268782993929651256519107625),
         ScalarType(-0.5384693101056830910363144207002088049672866069056),
         ScalarType(0.),
         ScalarType(0.5384693101056830910363144207002088049672866069056),
         ScalarType(0.9061798459386639927976268782993929651256519107625)});
  } else if constexpr (Order == 10) {
    return std::array<ScalarType, Order>(
        {ScalarType(-0.9739065285171717200779640120844520534282699466924),
         ScalarType(-0.8650633666889845107320966884234930485275430149653),
         ScalarType(-0.6794095682990244062343273651148735757692947118348),
         ScalarType(-0.4333953941292471907992659431657841622000718376562),
         ScalarType(-0.1488743389816312108848260011297199846175648594207),
         ScalarType(0.1488743389816312108848260011297199846175648594207),
         ScalarType(0.4333953941292471907992659431657841622000718376562),
         ScalarType(0.6794095682990244062343273651148735757692947118348),
         ScalarType(0.8650633666889845107320966884234930485275430149653),
         ScalarType(0.9739065285171717200779640120844520534282699466924)});

  } else if constexpr (Order == 50) {
    return std::array<ScalarType, Order>(
        {ScalarType(-0.9988664044200710501854594449742185059962435129041),
         ScalarType(-0.9940319694320907125851082004206947281574779710683),
         ScalarType(-0.9853540840480058823090096256324894040155926309454),
         ScalarType(-0.9728643851066920737133441046062520536691734070500),
         ScalarType(-0.9566109552428079429977456441566220940514341246260),
         ScalarType(-0.9366566189448779337808749472724966021537315980952),
         ScalarType(-0.9130785566557918930897356427716570947841881916978),
         ScalarType(-0.8859679795236130486375409824667536341942903107558),
         ScalarType(-0.8554297694299460846113626439347574676548330394869),
         ScalarType(-0.8215820708593359483562541108739395377607413833435),
         ScalarType(-0.7845558329003992639053051963409912008473162725564),
         ScalarType(-0.7444943022260685382605362526821942428701879313296),
         ScalarType(-0.7015524687068222510895462578836557281497192284847),
         ScalarType(-0.6558964656854393607816248640036798190414105283823),
         ScalarType(-0.6077029271849502391803817963918328936042050206760),
         ScalarType(-0.5571583045146500543155229096258016078158983822378),
         ScalarType(-0.5044581449074642016514591318491411926353786782707),
         ScalarType(-0.4498063349740387891471314677783758173150645134652),
         ScalarType(-0.3934143118975651273942292538238172702461394686727),
         ScalarType(-0.3355002454194373568369882572910716978412185930368),
         ScalarType(-0.2762881937795319903276452785211301857148015871320),
         ScalarType(-0.2160072368760417568472845326171013337057559725101),
         ScalarType(-0.1548905899981459020716286209411095012018502200549),
         ScalarType(-0.0931747015600861408544503776396003478856713839221),
         ScalarType(-0.0310983383271888761123289896659491942472962229600),
         ScalarType(0.0310983383271888761123289896659491942472962229600),
         ScalarType(0.0931747015600861408544503776396003478856713839221),
         ScalarType(0.1548905899981459020716286209411095012018502200549),
         ScalarType(0.2160072368760417568472845326171013337057559725101),
         ScalarType(0.2762881937795319903276452785211301857148015871320),
         ScalarType(0.3355002454194373568369882572910716978412185930368),
         ScalarType(0.3934143118975651273942292538238172702461394686727),
         ScalarType(0.4498063349740387891471314677783758173150645134652),
         ScalarType(0.5044581449074642016514591318491411926353786782707),
         ScalarType(0.5571583045146500543155229096258016078158983822378),
         ScalarType(0.6077029271849502391803817963918328936042050206760),
         ScalarType(0.6558964656854393607816248640036798190414105283823),
         ScalarType(0.7015524687068222510895462578836557281497192284847),
         ScalarType(0.7444943022260685382605362526821942428701879313296),
         ScalarType(0.7845558329003992639053051963409912008473162725564),
         ScalarType(0.8215820708593359483562541108739395377607413833435),
         ScalarType(0.8554297694299460846113626439347574676548330394869),
         ScalarType(0.8859679795236130486375409824667536341942903107558),
         ScalarType(0.9130785566557918930897356427716570947841881916978),
         ScalarType(0.9366566189448779337808749472724966021537315980952),
         ScalarType(0.9566109552428079429977456441566220940514341246260),
         ScalarType(0.9728643851066920737133441046062520536691734070500),
         ScalarType(0.9853540840480058823090096256324894040155926309454),
         ScalarType(0.9940319694320907125851082004206947281574779710683),
         ScalarType(0.9988664044200710501854594449742185059962435129041)});
  } else if constexpr (Order == 100) {
    return std::array<ScalarType, Order>(
        {ScalarType(-0.9997137267734412336782284693423006767183495273084),
         ScalarType(-0.9984919506395958184001633591863491623048548504206),
         ScalarType(-0.9962951347331251491861317322411310354364312881404),
         ScalarType(-0.9931249370374434596520098928487834707317714588665),
         ScalarType(-0.9889843952429917480044187458077366318393336371069),
         ScalarType(-0.9838775407060570154961001555110081673443670168508),
         ScalarType(-0.9778093584869182885537810884292019286352344942663),
         ScalarType(-0.9707857757637063319308978578975053885505571994782),
         ScalarType(-0.9628136542558155272936593260301663864373315067304),
         ScalarType(-0.9539007829254917428493369308943576446452214510147),
         ScalarType(-0.9440558701362559779627747064152187467397203733821),
         ScalarType(-0.9332885350430795459243336681308625040835460742970),
         ScalarType(-0.9216092981453339526669513284819874591245827977322),
         ScalarType(-0.9090295709825296904671263377891460644432772895846),
         ScalarType(-0.8955616449707269866985210224302277698481817689977),
         ScalarType(-0.8812186793850184155733168254278055824454944110213),
         ScalarType(-0.8660146884971646234107399696762429663803148393056),
         ScalarType(-0.8499645278795912842933625914201046540737907795007),
         ScalarType(-0.8330838798884008235429158338447556799074948303100),
         ScalarType(-0.8153892383391762543939887586492580053825503765137),
         ScalarType(-0.7968978923903144763895728821832459828895268559643),
         ScalarType(-0.7776279096494954756275513868344901065385397999361),
         ScalarType(-0.7575981185197071760356679644384007723131089719000),
         ScalarType(-0.7368280898020207055124277148201010028432784462471),
         ScalarType(-0.7153381175730564464599671227043659640843978385956),
         ScalarType(-0.6931491993558019659486479416754372655870000179307),
         ScalarType(-0.6702830156031410158025870143232266136698056840288),
         ScalarType(-0.6467619085141292798326303044586304350197337842485),
         ScalarType(-0.6226088602037077716041908451723122446538177322898),
         ScalarType(-0.5978474702471787212648065451493406363948991923205),
         ScalarType(-0.5725019326213811913168704435257254489600339496756),
         ScalarType(-0.5465970120650941674679942571817499039562417759375),
         ScalarType(-0.5201580198817630566468157494552085307689376904201),
         ScalarType(-0.4932107892081909335693087934493339909907233253586),
         ScalarType(-0.4657816497733580422492166233957545816116511102122),
         ScalarType(-0.4378974021720315131089780436221959621257017634841),
         ScalarType(-0.4095852916783015425288684000571577014953643891648),
         ScalarType(-0.3808729816246299567633625488695874037497072651237),
         ScalarType(-0.3517885263724217209723438295489705652493180963891),
         ScalarType(-0.3223603439005291517224765823983254274021916230231),
         ScalarType(-0.2926171880384719647375558882354943845615389891726),
         ScalarType(-0.2625881203715034791689293362549821411320226945355),
         ScalarType(-0.2323024818449739696495099632079641106975097715071),
         ScalarType(-0.2017898640957359972360488595303964629436920035590),
         ScalarType(-0.1710800805386032748875323747070898074658597251181),
         ScalarType(-0.1402031372361139732075146046824055166168730062634),
         ScalarType(-0.1091892035800611150034260065793848868848996299692),
         ScalarType(-0.0780685828134366366948173712015525739763500274485),
         ScalarType(-0.0468716824215916316149239129338483095370653990860),
         ScalarType(-0.0156289844215430828722166999974293401477561828556),
         ScalarType(0.0156289844215430828722166999974293401477561828556),
         ScalarType(0.0468716824215916316149239129338483095370653990860),
         ScalarType(0.0780685828134366366948173712015525739763500274485),
         ScalarType(0.1091892035800611150034260065793848868848996299692),
         ScalarType(0.1402031372361139732075146046824055166168730062634),
         ScalarType(0.1710800805386032748875323747070898074658597251181),
         ScalarType(0.2017898640957359972360488595303964629436920035590),
         ScalarType(0.2323024818449739696495099632079641106975097715071),
         ScalarType(0.2625881203715034791689293362549821411320226945355),
         ScalarType(0.2926171880384719647375558882354943845615389891726),
         ScalarType(0.3223603439005291517224765823983254274021916230231),
         ScalarType(0.3517885263724217209723438295489705652493180963891),
         ScalarType(0.3808729816246299567633625488695874037497072651237),
         ScalarType(0.4095852916783015425288684000571577014953643891648),
         ScalarType(0.4378974021720315131089780436221959621257017634841),
         ScalarType(0.4657816497733580422492166233957545816116511102122),
         ScalarType(0.4932107892081909335693087934493339909907233253586),
         ScalarType(0.5201580198817630566468157494552085307689376904201),
         ScalarType(0.5465970120650941674679942571817499039562417759375),
         ScalarType(0.5725019326213811913168704435257254489600339496756),
         ScalarType(0.5978474702471787212648065451493406363948991923205),
         ScalarType(0.6226088602037077716041908451723122446538177322898),
         ScalarType(0.6467619085141292798326303044586304350197337842485),
         ScalarType(0.6702830156031410158025870143232266136698056840288),
         ScalarType(0.6931491993558019659486479416754372655870000179307),
         ScalarType(0.7153381175730564464599671227043659640843978385956),
         ScalarType(0.7368280898020207055124277148201010028432784462471),
         ScalarType(0.7575981185197071760356679644384007723131089719000),
         ScalarType(0.7776279096494954756275513868344901065385397999361),
         ScalarType(0.7968978923903144763895728821832459828895268559643),
         ScalarType(0.8153892383391762543939887586492580053825503765137),
         ScalarType(0.8330838798884008235429158338447556799074948303100),
         ScalarType(0.8499645278795912842933625914201046540737907795007),
         ScalarType(0.8660146884971646234107399696762429663803148393056),
         ScalarType(0.8812186793850184155733168254278055824454944110213),
         ScalarType(0.8955616449707269866985210224302277698481817689977),
         ScalarType(0.9090295709825296904671263377891460644432772895846),
         ScalarType(0.9216092981453339526669513284819874591245827977322),
         ScalarType(0.9332885350430795459243336681308625040835460742970),
         ScalarType(0.9440558701362559779627747064152187467397203733821),
         ScalarType(0.9539007829254917428493369308943576446452214510147),
         ScalarType(0.9628136542558155272936593260301663864373315067304),
         ScalarType(0.9707857757637063319308978578975053885505571994782),
         ScalarType(0.9778093584869182885537810884292019286352344942663),
         ScalarType(0.9838775407060570154961001555110081673443670168508),
         ScalarType(0.9889843952429917480044187458077366318393336371069),
         ScalarType(0.9931249370374434596520098928487834707317714588665),
         ScalarType(0.9962951347331251491861317322411310354364312881404),
         ScalarType(0.9984919506395958184001633591863491623048548504206),
         ScalarType(0.9997137267734412336782284693423006767183495273084)});
  }
}

template <class ScalarType, UnsignedIndex_t Order>
const std::array<ScalarType, Order> WeightsGauss(void) {
  if constexpr (Order == 5) {
    return std::array<ScalarType, Order>(
        {ScalarType(0.23692688505618908751426404071991736264326000221241),
         ScalarType(0.47862867049936646804129151483563819291229555334314),
         ScalarType(0.5688888888888888888888888888888888888888888888889),
         ScalarType(0.47862867049936646804129151483563819291229555334314),
         ScalarType(0.23692688505618908751426404071991736264326000221241)});
  } else if constexpr (Order == 10) {
    return std::array<ScalarType, Order>(
        {ScalarType(0.06667134430868813759356880989333179285786483432016),
         ScalarType(0.14945134915058059314577633965769733240255663966943),
         ScalarType(0.21908636251598204399553493422816319245877187052268),
         ScalarType(0.26926671930999635509122692156946935285975993846088),
         ScalarType(0.29552422471475287017389299465133832942104671702685),
         ScalarType(0.29552422471475287017389299465133832942104671702685),
         ScalarType(0.26926671930999635509122692156946935285975993846088),
         ScalarType(0.21908636251598204399553493422816319245877187052268),
         ScalarType(0.14945134915058059314577633965769733240255663966943),
         ScalarType(0.06667134430868813759356880989333179285786483432016)});
  } else if constexpr (Order == 50) {
    return std::array<ScalarType, Order>(
        {ScalarType(0.0029086225531551409584007243428554808066729964599463),
         ScalarType(0.006759799195745401502778878177985031801873832406467),
         ScalarType(0.010590548383650969263569681499241022339401819086459),
         ScalarType(0.014380822761485574419378908927324349937031786170588),
         ScalarType(0.018115560713489390351259943422354619844667317049734),
         ScalarType(0.021780243170124792981592069062690341227313462357934),
         ScalarType(0.025360673570012390440194878385442723460161259975712),
         ScalarType(0.028842993580535198029906373113232432517846865593537),
         ScalarType(0.032213728223578016648165827323003953448589058833425),
         ScalarType(0.035459835615146154160734611000975797096960000496984),
         ScalarType(0.038568756612587675244770150236385934864771705000519),
         ScalarType(0.041528463090147697422411978964067017808977975485841),
         ScalarType(0.044327504338803275492022286830394197460761298355453),
         ScalarType(0.046955051303948432965633013634987682514064306186053),
         ScalarType(0.049400938449466314921243580751432728692287050966613),
         ScalarType(0.05165570306958113848990529584009527964982544939544),
         ScalarType(0.05371062188899624652345879725566455276802321352992),
         ScalarType(0.05555774480621251762356742561226949759513529998390),
         ScalarType(0.05718992564772838372302931506599316301157537225709),
         ScalarType(0.05860084981322244583512243663084846620976751344403),
         ScalarType(0.05978505870426545750957640531258523079666604207267),
         ScalarType(0.06073797084177021603175001538481100160979927323540),
         ScalarType(0.06145589959031666375640678608391537509726757576401),
         ScalarType(0.06193606742068324338408750978083068857287705669124),
         ScalarType(0.06217661665534726232103310736061343086768246920103),
         ScalarType(0.06217661665534726232103310736061343086768246920103),
         ScalarType(0.06193606742068324338408750978083068857287705669124),
         ScalarType(0.06145589959031666375640678608391537509726757576401),
         ScalarType(0.06073797084177021603175001538481100160979927323540),
         ScalarType(0.05978505870426545750957640531258523079666604207267),
         ScalarType(0.05860084981322244583512243663084846620976751344403),
         ScalarType(0.05718992564772838372302931506599316301157537225709),
         ScalarType(0.05555774480621251762356742561226949759513529998390),
         ScalarType(0.05371062188899624652345879725566455276802321352992),
         ScalarType(0.05165570306958113848990529584009527964982544939544),
         ScalarType(0.049400938449466314921243580751432728692287050966613),
         ScalarType(0.046955051303948432965633013634987682514064306186053),
         ScalarType(0.044327504338803275492022286830394197460761298355453),
         ScalarType(0.041528463090147697422411978964067017808977975485841),
         ScalarType(0.038568756612587675244770150236385934864771705000519),
         ScalarType(0.035459835615146154160734611000975797096960000496984),
         ScalarType(0.032213728223578016648165827323003953448589058833425),
         ScalarType(0.028842993580535198029906373113232432517846865593537),
         ScalarType(0.025360673570012390440194878385442723460161259975712),
         ScalarType(0.021780243170124792981592069062690341227313462357934),
         ScalarType(0.018115560713489390351259943422354619844667317049734),
         ScalarType(0.014380822761485574419378908927324349937031786170588),
         ScalarType(0.010590548383650969263569681499241022339401819086459),
         ScalarType(0.006759799195745401502778878177985031801873832406467),
         ScalarType(0.0029086225531551409584007243428554808066729964599463)});
  } else if constexpr (Order == 100) {
    return std::array<ScalarType, Order>(
        {ScalarType(0.0007346344905056717304063206583303363906704735624829),
         ScalarType(0.0017093926535181052395293583714911952437313854914626),
         ScalarType(0.0026839253715534824194395904290011200819311149509983),
         ScalarType(0.0036559612013263751823424587275251956992065674051522),
         ScalarType(0.0046244500634221193510957890829784766503524952948945),
         ScalarType(0.005588428003865515157211946348439210731318694008077),
         ScalarType(0.006546948450845322764152103331495263699938363366476),
         ScalarType(0.007499073255464711578828744016397783163583478948145),
         ScalarType(0.008443871469668971402620834902301001934644459884101),
         ScalarType(0.009380419653694457951418237660812118730787043238674),
         ScalarType(0.010307802574868969585782101727835377976058343841426),
         ScalarType(0.011225114023185977117221573366333584777226419564382),
         ScalarType(0.012131457662979497407744792448748170736963123311261),
         ScalarType(0.013025947892971542285558583758901790134964735841750),
         ScalarType(0.013907710703718772687954149108004637795180812143119),
         ScalarType(0.014775884527441301768879987520354257169388743114604),
         ScalarType(0.015629621077546002723936865953791925552469979809937),
         ScalarType(0.016468086176145212643104980088210780821167661603798),
         ScalarType(0.017290460568323582439344198366741674811623508565168),
         ScalarType(0.018095940722128116664390751420493031347578744958390),
         ScalarType(0.018883739613374904552941165881543234297111276347424),
         ScalarType(0.019653087494435305865381470245444065555269599491251),
         ScalarType(0.020403232646209432766838851657583770605709699302620),
         ScalarType(0.021133442112527641542672300440969681635329728874515),
         ScalarType(0.021843002416247386313953741304398024765348999823246),
         ScalarType(0.022531220256336272701796970931673962340158935348709),
         ScalarType(0.023197423185254121622488854182727288451154485736087),
         ScalarType(0.023840960265968205962560411902283432144707449092617),
         ScalarType(0.024461202707957052719975023349772890646295732397804),
         ScalarType(0.025057544481579589703764225620923264223838558527929),
         ScalarType(0.025629402910208116075642009862150870926977670020268),
         ScalarType(0.026176219239545676342308741757301885011275131190695),
         ScalarType(0.026697459183570962660384664186336350634655750390012),
         ScalarType(0.027192613446576880136491567802170692266987896012001),
         ScalarType(0.027661198220792388294204155870426455292400358664220),
         ScalarType(0.028102755659101173317648330186994550451418099400205),
         ScalarType(0.028516854322395097990936762864457873259842725483969),
         ScalarType(0.028903089601125203134876228134515265315607868055261),
         ScalarType(0.029261084110638276620119023495640954443084195045348),
         ScalarType(0.029590488059912642511754510678836585172806285071368),
         ScalarType(0.029890979593332830916836806668595827658091414260799),
         ScalarType(0.030162265105169144919068681610479232657102325782706),
         ScalarType(0.030404079526454820016507859818825176605607248310117),
         ScalarType(0.030616186583980448496459443262053192853086023789058),
         ScalarType(0.030798379031152590427713903030559760094970834470371),
         ScalarType(0.030950478850490988234063463470747927382987177766941),
         ScalarType(0.031072337427566516587810170242918034845915436347957),
         ScalarType(0.031163835696209906783818321217186653343836368683928),
         ScalarType(0.031224884254849357732376498648098134881802740682184),
         ScalarType(0.031255423453863356947642474386198028787833836726091),
         ScalarType(0.031255423453863356947642474386198028787833836726091),
         ScalarType(0.031224884254849357732376498648098134881802740682184),
         ScalarType(0.031163835696209906783818321217186653343836368683928),
         ScalarType(0.031072337427566516587810170242918034845915436347957),
         ScalarType(0.030950478850490988234063463470747927382987177766941),
         ScalarType(0.030798379031152590427713903030559760094970834470371),
         ScalarType(0.030616186583980448496459443262053192853086023789058),
         ScalarType(0.030404079526454820016507859818825176605607248310117),
         ScalarType(0.030162265105169144919068681610479232657102325782706),
         ScalarType(0.029890979593332830916836806668595827658091414260799),
         ScalarType(0.029590488059912642511754510678836585172806285071368),
         ScalarType(0.029261084110638276620119023495640954443084195045348),
         ScalarType(0.028903089601125203134876228134515265315607868055261),
         ScalarType(0.028516854322395097990936762864457873259842725483969),
         ScalarType(0.028102755659101173317648330186994550451418099400205),
         ScalarType(0.027661198220792388294204155870426455292400358664220),
         ScalarType(0.027192613446576880136491567802170692266987896012001),
         ScalarType(0.026697459183570962660384664186336350634655750390012),
         ScalarType(0.026176219239545676342308741757301885011275131190695),
         ScalarType(0.025629402910208116075642009862150870926977670020268),
         ScalarType(0.025057544481579589703764225620923264223838558527929),
         ScalarType(0.024461202707957052719975023349772890646295732397804),
         ScalarType(0.023840960265968205962560411902283432144707449092617),
         ScalarType(0.023197423185254121622488854182727288451154485736087),
         ScalarType(0.022531220256336272701796970931673962340158935348709),
         ScalarType(0.021843002416247386313953741304398024765348999823246),
         ScalarType(0.021133442112527641542672300440969681635329728874515),
         ScalarType(0.020403232646209432766838851657583770605709699302620),
         ScalarType(0.019653087494435305865381470245444065555269599491251),
         ScalarType(0.018883739613374904552941165881543234297111276347424),
         ScalarType(0.018095940722128116664390751420493031347578744958390),
         ScalarType(0.017290460568323582439344198366741674811623508565168),
         ScalarType(0.016468086176145212643104980088210780821167661603798),
         ScalarType(0.015629621077546002723936865953791925552469979809937),
         ScalarType(0.014775884527441301768879987520354257169388743114604),
         ScalarType(0.013907710703718772687954149108004637795180812143119),
         ScalarType(0.013025947892971542285558583758901790134964735841750),
         ScalarType(0.012131457662979497407744792448748170736963123311261),
         ScalarType(0.011225114023185977117221573366333584777226419564382),
         ScalarType(0.010307802574868969585782101727835377976058343841426),
         ScalarType(0.009380419653694457951418237660812118730787043238674),
         ScalarType(0.008443871469668971402620834902301001934644459884101),
         ScalarType(0.007499073255464711578828744016397783163583478948145),
         ScalarType(0.006546948450845322764152103331495263699938363366476),
         ScalarType(0.005588428003865515157211946348439210731318694008077),
         ScalarType(0.0046244500634221193510957890829784766503524952948945),
         ScalarType(0.0036559612013263751823424587275251956992065674051522),
         ScalarType(0.0026839253715534824194395904290011200819311149509983),
         ScalarType(0.0017093926535181052395293583714911952437313854914626),
         ScalarType(0.0007346344905056717304063206583303363906704735624829)});
  }
}

/* This returns the algebraic signed distance to a paraboloid
 * returns < 0 number if a_pt is below the paraboloid
 * returns > 0 number if a_pt is above the paraboloid */
template <class ScalarType>
inline ScalarType signedDistance(
    const PtBase<ScalarType>& a_pt,
    const AlignedParaboloidBase<ScalarType>& a_paraboloid) {
  return a_paraboloid.a() * a_pt[0] * a_pt[0] +
         a_paraboloid.b() * a_pt[1] * a_pt[1] + a_pt[2];
}

template <class ScalarType>
inline ScalarType CstHalf(void) {
  return ScalarType(0.5);
}

template <class ScalarType>
inline ScalarType CstThird(void) {
  return ScalarType(0.333333333333333333333333333333);
}

template <class ScalarType>
inline ScalarType CstFourth(void) {
  return ScalarType(0.25);
}

template <class ScalarType>
inline ScalarType CstFifth(void) {
  return ScalarType(0.2);
}

template <class ScalarType>
inline ScalarType CstSixth(void) {
  return ScalarType(0.166666666666666666666666666667);
}

template <class ScalarType>
inline ScalarType CstSeventh(void) {
  return ScalarType(1) / ScalarType(7);
}

template <class ScalarType>
inline ScalarType CstEighth(void) {
  return ScalarType(0.125);
}
template <class ScalarType, UnsignedIndex_t OrderX, UnsignedIndex_t OrderY,
          UnsignedIndex_t OrderZ, UnsignedIndex_t ProjDir>
inline ScalarType MomentPlaneIntegrand(const PtBase<ScalarType>& a_position,
                                       const PtBase<ScalarType>& a_derivative,
                                       const NormalBase<ScalarType>& a_normal,
                                       const ScalarType& d,
                                       const ScalarType& w) {
  const ScalarType &x = a_position[0], &y = a_position[1], &z = a_position[2];
  const ScalarType &dy = a_derivative[1], &dz = a_derivative[2];
  const ScalarType &nx = a_normal[0], &ny = a_normal[1], &nz = a_normal[2];
  const ScalarType inv2 = CstHalf<ScalarType>();
  const ScalarType inv3 = CstThird<ScalarType>();
  const ScalarType inv4 = CstFourth<ScalarType>();
  const ScalarType inv6 = CstSixth<ScalarType>();
  const ScalarType inv8 = CstEighth<ScalarType>();
  // M0
  if constexpr (OrderX == 0 && OrderY == 0 && OrderZ == 0) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz * y * z;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz * x * z;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      return w * dy * (d * x - inv2 * nx * x * x - ny * x * y);
    }
  }
  // M1x
  else if constexpr (OrderX == 1 && OrderY == 0 && OrderZ == 0) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz * (d * y * z - inv2 * ny * y * y * z - nz * y * z * z) / nx;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz * inv2 * x * x * z;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      const auto x2 = x * x;
      return w * dy * (inv2 * d * x2 - inv3 * nx * x * x2 - inv2 * ny * x2 * y);
    }
  }
  // M1y
  else if constexpr (OrderX == 0 && OrderY == 1 && OrderZ == 0) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz * inv2 * y * y * z;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz * (d * x * z - inv2 * nx * x * x * z - nz * x * z * z) / ny;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      return w * dy * (d * x * y - inv2 * nx * x * x * y - ny * x * y * y);
    }
  }
  // M1z
  else if constexpr (OrderX == 0 && OrderY == 0 && OrderZ == 1) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz * inv2 * y * z * z;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz * inv2 * x * z * z;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      const auto x2 = x * x;
      return w * dy *
             (inv2 * d * d * x - inv2 * d * nx * x2 + inv6 * nx * nx * x * x2 -
              d * ny * x * y + inv2 * nx * ny * x2 * y +
              inv2 * ny * ny * x * y * y) /
             nz;
    }
  }
  // M2xx
  else if constexpr (OrderX == 2 && OrderY == 0 && OrderZ == 0) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz *
             (d * d * y * z - d * ny * y * y * z +
              inv3 * ny * ny * y * y * y * z -
              ScalarType(2) * d * nz * y * z * z + ny * nz * y * y * z * z +
              nz * nz * y * z * z * z) /
             (nx * nx);
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz * inv3 * x * x * x * z;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      const auto x3 = x * x * x;
      return w * dy * (inv3 * d * x3 - inv4 * nx * x * x3 - inv3 * ny * x3 * y);
    }
  }
  // M2xy
  else if constexpr (OrderX == 1 && OrderY == 1 && OrderZ == 0) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz *
             (inv2 * d * y * y * z - inv3 * ny * y * y * y * z -
              inv2 * nz * y * y * z * z) /
             nx;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz *
             (inv2 * d * x * x * z - inv3 * nx * x * x * x * z -
              inv2 * nz * x * x * z * z) /
             ny;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      const auto x2 = x * x;
      return w * dy *
             (inv2 * d * x2 * y - inv3 * nx * x * x2 * y -
              inv2 * ny * x2 * y * y);
    }
  }
  // M2xz
  else if constexpr (OrderX == 1 && OrderY == 0 && OrderZ == 1) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz *
             (inv2 * d * y * z * z - inv4 * ny * y * y * z * z -
              inv2 * nz * y * z * z * z) /
             nx;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz * inv4 * x * x * z * z;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      const auto x2 = x * x;
      const auto x3 = x * x2;
      return w * dy *
             (inv4 * d * d * x2 - inv3 * d * nx * x3 + inv8 * nx * nx * x * x3 -
              inv2 * d * ny * x2 * y + inv3 * nx * ny * x3 * y +
              inv4 * ny * ny * x2 * y * y) /
             nz;
    }
  }
  // M2yy
  else if constexpr (OrderX == 0 && OrderY == 2 && OrderZ == 0) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz * inv3 * y * y * y * z;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz *
             (d * d * x * z - d * nx * x * x * z +
              inv3 * nx * nx * x * x * x * z -
              ScalarType(2) * d * nz * x * z * z + nx * nz * x * x * z * z +
              nz * nz * x * z * z * z) /
             (ny * ny);
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      const auto y2 = y * y;
      return w * dy * (d * x * y2 - inv2 * nx * x * x * y2 - ny * x * y * y2);
    }
  }
  // M2yz
  else if constexpr (OrderX == 0 && OrderY == 1 && OrderZ == 1) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz * inv4 * y * y * z * z;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz *
             (inv2 * d * x * z * z - inv4 * nx * x * x * z * z -
              inv2 * nz * x * z * z * z) /
             ny;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      const auto x2 = x * x;
      const auto y2 = y * y;
      return w * dy *
             (inv2 * d * d * x * y - inv2 * d * nx * x2 * y +
              inv6 * nx * nx * x2 * x * y - d * ny * x * y2 +
              inv2 * nx * ny * x2 * y2 + inv2 * ny * ny * x * y * y2) /
             nz;
    }
  }
  // M2zz
  else if constexpr (OrderX == 0 && OrderY == 0 && OrderZ == 2) {
    // Projection on YZ
    if constexpr (ProjDir == 0) {
      return w * dz * inv3 * y * z * z * z;
    }
    // Projection on ZX
    else if constexpr (ProjDir == 1) {
      return w * dz * inv3 * x * z * z * z;
    }
    // Projection on XY
    else if constexpr (ProjDir == 2) {
      const auto d2 = d * d;
      const auto x2 = x * x;
      const auto x3 = x * x2;
      const auto y2 = y * y;
      const auto nx2 = nx * nx;
      const auto ny2 = ny * ny;
      return w * dy *
             (inv3 * d2 * d * x - inv2 * d2 * nx * x2 + inv3 * d * nx2 * x3 -
              inv2 * inv6 * nx * nx2 * x * x3 - d2 * ny * x * y +
              d * nx * ny * x2 * y - inv3 * nx2 * ny * x3 * y +
              d * ny2 * x * y2 - inv2 * nx * ny2 * x2 * y2 -
              inv3 * ny * ny2 * x * y * y2) /
             (nz * nz);
    }
  } else {
    return ScalarType(0);
  }
}

template <class ScalarType, UnsignedIndex_t OrderX, UnsignedIndex_t OrderY,
          UnsignedIndex_t OrderZ>
inline ScalarType MomentParaboloidIntegrand(
    const PtBase<ScalarType>& a_position,
    const PtBase<ScalarType>& a_derivative, const ScalarType& A,
    const ScalarType& B) {
  const ScalarType &x = a_position[0], &y = a_position[1],
                   &dy = a_derivative[1];
  const ScalarType inv2 = CstHalf<ScalarType>();
  const ScalarType inv3 = CstThird<ScalarType>();
  const ScalarType inv4 = CstFourth<ScalarType>();
  const ScalarType inv5 = CstFifth<ScalarType>();
  const ScalarType inv6 = CstSixth<ScalarType>();
  const ScalarType inv7 = CstSeventh<ScalarType>();

  // M0
  if constexpr (OrderX == 0 && OrderY == 0 && OrderZ == 0) {
    return -dy * (inv3 * A * x * x * x + B * x * y * y);
  }
  // M1x
  else if constexpr (OrderX == 1 && OrderY == 0 && OrderZ == 0) {
    const auto x2 = x * x;
    return -dy * (inv4 * A * x2 * x2 + inv2 * B * x2 * y * y);
  }
  // M1y
  else if constexpr (OrderX == 0 && OrderY == 1 && OrderZ == 0) {
    return -dy * (inv3 * A * x * x * x * y + B * x * y * y * y);
  }
  // M1z
  else if constexpr (OrderX == 0 && OrderY == 0 && OrderZ == 1) {
    const auto y2 = y * y;
    const auto x3 = x * x * x;
    return dy * (inv2 * inv5 * A * A * x3 * x * x + inv3 * A * B * x3 * y2 +
                 inv2 * B * B * x * y2 * y2);
  }
  // M2xx
  else if constexpr (OrderX == 2 && OrderY == 0 && OrderZ == 0) {
    const auto x3 = x * x * x;
    return -dy * (inv5 * A * x3 * x * x + inv3 * B * x3 * y * y);
  }
  // M2xy
  else if constexpr (OrderX == 1 && OrderY == 1 && OrderZ == 0) {
    const auto x2 = x * x;
    return -dy * (inv4 * A * x2 * x2 * y + inv2 * B * x2 * y * y * y);
  }
  // M2xz
  else if constexpr (OrderX == 1 && OrderY == 0 && OrderZ == 1) {
    const auto x2 = x * x;
    const auto y2 = y * y;
    const auto x4 = x2 * x2;
    return dy * (inv4 * inv3 * A * A * x4 * x2 * y + inv4 * A * B * x4 * y2 +
                 inv4 * B * B * x2 * y2 * y2);
  }
  // M2yy
  else if constexpr (OrderX == 0 && OrderY == 2 && OrderZ == 0) {
    const auto y2 = y * y;
    return -dy * (inv3 * A * x * x * x * y2 + B * x * y2 * y2);
  }
  // M2yz
  else if constexpr (OrderX == 0 && OrderY == 1 && OrderZ == 1) {
    const auto x3 = x * x * x;
    const auto y3 = y * y * y;
    return dy * (inv5 * inv2 * A * A * x3 * x * x * y + inv3 * A * B * x3 * y3 +
                 inv2 * B * B * x * y3 * y * y);
  }
  // M2zz
  else if constexpr (OrderX == 0 && OrderY == 0 && OrderZ == 2) {
    const auto x2 = x * x;
    const auto x3 = x * x2;
    const auto x5 = x2 * x3;
    const auto y2 = y * y;
    const auto y4 = y2 * y2;
    return -dy *
           (inv3 * inv7 * A * A * A * x5 * x2 + inv5 * A * A * B * x5 * y2 +
            inv3 * A * B * B * x3 * y4 + inv3 * B * B * B * x * y4 * y2);
  }
  return ScalarType(0);
}

template <UnsignedIndex_t ProjDir, class ReturnType, class ScalarType,
          class MomentFunctorType>
inline ReturnType MomentsIntegrandArc(const ScalarType a_t,
                                      const MomentFunctorType& a_functor) {
  using ReturnScalarType = typename ReturnType::value_type;
  const auto der_t = a_functor.der_t(a_t);
  const auto pos_t = a_functor.pos_t(a_t);
  const auto& normal = a_functor.face_normal();
  const auto& dist = a_functor.face_distance();
  const auto& weight = a_functor.norm_weight();
  const auto& a = a_functor.a();
  const auto& b = a_functor.b();
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 0, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight)));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    auto integrand = ReturnType::fromScalarConstant(ReturnScalarType(0));
    integrand.volume() = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 0, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[0] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 1, 0, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 1, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[1] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 1, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 1, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand.centroid()[2] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 0, 1>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 1, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    return integrand;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    auto integrand = ReturnType::fromScalarConstant(ReturnScalarType(0));
    integrand[0] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 0, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[1] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 1, 0, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 1, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[2] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 1, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 1, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[3] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 0, 1>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 1, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[4] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 2, 0, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 2, 0, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[5] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 1, 1, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 1, 1, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[6] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 1, 0, 1>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 1, 0, 1, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[7] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 2, 0>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 2, 0, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[8] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 1, 1>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 1, 1, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    integrand[9] = ReturnScalarType(
        MomentParaboloidIntegrand<ScalarType, 0, 0, 2>(pos_t, der_t, a, b) -
        MomentPlaneIntegrand<ScalarType, 0, 0, 2, ProjDir>(pos_t, der_t, normal,
                                                           dist, weight));
    return integrand;
  } else {
    std::cout << "Paraboloid M>2 not available yet" << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

template <UnsignedIndex_t ProjDir, class ReturnType, class ScalarType,
          class MomentFunctorType>
inline ReturnType MomentsIntegrandLine(const ScalarType a_t,
                                       const MomentFunctorType& a_functor) {
  using ReturnScalarType = typename ReturnType::value_type;
  const auto der_t = a_functor.der_t(a_t);
  const auto pos_t = a_functor.pos_t(a_t);
  const auto& normal = a_functor.face_normal();
  const auto& dist = a_functor.face_distance();
  const auto& weight = a_functor.norm_weight();
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight)));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    auto integrand = ReturnType::fromScalarConstant(ReturnScalarType(0));
    integrand.volume() =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand.centroid()[0] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 1, 0, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand.centroid()[1] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 1, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand.centroid()[2] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 0, 1, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    return integrand;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    auto integrand = ReturnType::fromScalarConstant(ReturnScalarType(0));
    integrand[0] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 0, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[1] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 1, 0, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[2] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 1, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[3] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 0, 1, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[4] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 2, 0, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[5] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 1, 1, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[6] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 1, 0, 1, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[7] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 2, 0, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[8] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 1, 1, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    integrand[9] =
        ReturnScalarType(MomentPlaneIntegrand<ScalarType, 0, 0, 2, ProjDir>(
            pos_t, der_t, normal, dist, weight));
    return integrand;
  } else {
    std::cout << "Paraboloid M>2 not available yet" << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

template <class ReturnType, class ScalarType, UnsignedIndex_t QuadRuleOrder>
class ParaboloidMomentLineIntegrator {
 public:
  using ReturnScalarType = typename ReturnType::value_type;
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  ParaboloidMomentLineIntegrator(const Pt& a_pt_start, const Pt& a_pt_end,
                                 const Normal& a_face_normal,
                                 const UnsignedIndex_t a_proj_dir)
      : lower_limit_m(ScalarType(0)),
        upper_limit_m(ScalarType(1)),
        start_m(a_pt_start),
        end_m(a_pt_end),
        face_normal_m(a_face_normal),
        face_distance_m(a_face_normal * a_pt_start) {
    switch (a_proj_dir) {
      case 0:
        integrand_m = MomentsIntegrandLine<0>;
        norm_weight_m = face_normal_m[2] / face_normal_m[0];
        break;
      case 1:
        integrand_m = MomentsIntegrandLine<1>;
        norm_weight_m = -face_normal_m[2] / face_normal_m[1];
        break;
      case 2:
        integrand_m = MomentsIntegrandLine<2>;
        norm_weight_m = ScalarType(1) / face_normal_m[2];
        break;
      default:
        UnsignedIndex_t max_component_index = 0;
        ScalarType max_component = fabs(face_normal_m[0]);
        for (UnsignedIndex_t d = 1; d < 3; ++d) {
          if (fabs(face_normal_m[d]) > max_component) {
            max_component_index = d;
            max_component = fabs(face_normal_m[d]);
          }
        }
        switch (max_component_index) {
          case 0:
            integrand_m = MomentsIntegrandLine<0>;
            norm_weight_m = face_normal_m[2] / face_normal_m[0];
            break;
          case 1:
            integrand_m = MomentsIntegrandLine<1>;
            norm_weight_m = -face_normal_m[2] / face_normal_m[1];
            break;
          default:
            integrand_m = MomentsIntegrandLine<2>;
            norm_weight_m = ScalarType(1) / face_normal_m[2];
        }
    }
  }
  const ReturnType operator()(const ScalarType a_t) const {
    return (*integrand_m)(a_t, (*this));
  }
  const ReturnType integrate(void) {
    const auto one = ScalarType(1), inv2 = ScalarType(0.5);
    const auto& abscissea = AbscissaeGauss<ScalarType, QuadRuleOrder>();
    const auto& weights = WeightsGauss<ScalarType, QuadRuleOrder>();
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(0));
    for (UnsignedIndex_t i = 0; i < QuadRuleOrder; ++i) {
      const auto t = lower_limit_m + (upper_limit_m - lower_limit_m) * inv2 *
                                         (one + abscissea[i]);
      moments += ReturnScalarType(inv2 * weights[i]) * (*this)(t);
    }
    return moments;
  }
  inline const Pt der_t(const ScalarType a_t) const {
    return Pt(end_m - start_m);
  }
  inline const Pt pos_t(const ScalarType a_t) const {
    return Pt(a_t * end_m + (ScalarType(1) - a_t) * start_m);
  }
  inline const Normal& face_normal(void) const { return face_normal_m; }
  inline const ScalarType& face_distance(void) const { return face_distance_m; }
  inline const ScalarType& norm_weight(void) const { return norm_weight_m; }

 private:
  const Pt start_m, end_m;
  const Normal face_normal_m;
  const ScalarType face_distance_m;
  ScalarType lower_limit_m, upper_limit_m;
  ScalarType norm_weight_m;
  ReturnType (*integrand_m)(const ScalarType,
                            const ParaboloidMomentLineIntegrator&);
};

template <class ReturnType, class ScalarType, UnsignedIndex_t QuadRuleOrder>
class ParaboloidMomentArcIntegrator {
 public:
  using ReturnScalarType = typename ReturnType::value_type;
  using AlignedParaboloid = AlignedParaboloidBase<ScalarType>;
  using RationalBezierArc = RationalBezierArcBase<ScalarType>;
  using Pt = PtBase<ScalarType>;
  using Normal = NormalBase<ScalarType>;
  ParaboloidMomentArcIntegrator(const AlignedParaboloid& a_paraboloid,
                                const RationalBezierArc& a_arc,
                                const Normal& a_face_normal,
                                const UnsignedIndex_t a_proj_dir)
      : lower_limit_m(ScalarType(0)),
        upper_limit_m(ScalarType(1)),
        arc_m(a_arc),
        face_normal_m(a_face_normal),
        a_m(a_paraboloid.a()),
        b_m(a_paraboloid.b()),
        face_distance_m(a_face_normal * a_arc.start_point()) {
    switch (a_proj_dir) {
      case 0:
        integrand_m = MomentsIntegrandArc<0>;
        norm_weight_m = face_normal_m[2] / face_normal_m[0];
        break;
      case 1:
        integrand_m = MomentsIntegrandArc<1>;
        norm_weight_m = -face_normal_m[2] / face_normal_m[1];
        break;
      case 2:
        integrand_m = MomentsIntegrandArc<2>;
        norm_weight_m = ScalarType(1) / face_normal_m[2];
        break;
      default:
        UnsignedIndex_t max_component_index = 0;
        ScalarType max_component = fabs(face_normal_m[0]);
        for (UnsignedIndex_t d = 1; d < 3; ++d) {
          if (fabs(face_normal_m[d]) > max_component) {
            max_component_index = d;
            max_component = fabs(face_normal_m[d]);
          }
        }
        switch (max_component_index) {
          case 0:
            integrand_m = MomentsIntegrandArc<0>;
            norm_weight_m = face_normal_m[2] / face_normal_m[0];
            break;
          case 1:
            integrand_m = MomentsIntegrandArc<1>;
            norm_weight_m = -face_normal_m[2] / face_normal_m[1];
            break;
          default:
            integrand_m = MomentsIntegrandArc<2>;
            norm_weight_m = ScalarType(1) / face_normal_m[2];
        }
    }
  }
  const ReturnType operator()(const ScalarType a_t) const {
    return (*integrand_m)(a_t, (*this));
  }
  const ReturnType integrate(void) {
    const auto one = ScalarType(1), inv2 = ScalarType(0.5);
    const auto& abscissea = AbscissaeGauss<ScalarType, QuadRuleOrder>();
    const auto& weights = WeightsGauss<ScalarType, QuadRuleOrder>();
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(0));
    for (UnsignedIndex_t i = 0; i < QuadRuleOrder; ++i) {
      const auto t = lower_limit_m + (upper_limit_m - lower_limit_m) * inv2 *
                                         (one + abscissea[i]);
      moments += ReturnScalarType(inv2 * weights[i]) * (*this)(t);
    }
    return moments;
  }
  inline const Pt der_t(const ScalarType a_t) const {
    return arc_m.derivative(a_t);
  }
  inline const Pt pos_t(const ScalarType a_t) const { return arc_m.point(a_t); }
  inline const Normal& face_normal(void) const { return face_normal_m; }
  inline const ScalarType& face_distance(void) const { return face_distance_m; }
  inline const ScalarType& norm_weight(void) const { return norm_weight_m; }
  inline const ScalarType& a(void) const { return a_m; }
  inline const ScalarType& b(void) const { return b_m; }

 private:
  const RationalBezierArc arc_m;
  const Normal face_normal_m;
  const ScalarType face_distance_m, a_m, b_m;
  ScalarType norm_weight_m;
  ScalarType lower_limit_m, upper_limit_m;
  ReturnType (*integrand_m)(const ScalarType,
                            const ParaboloidMomentArcIntegrator&);
};

/******************************************************************************/
/*********************** First moment contribution
 * ****************************/
/******************************************************************************/
/* This compute the first contribution to the moments (arising from the
 * integration of the face plane primitives on the poligonized clipped faces)
 */
template <class ReturnType, class ScalarType>
ReturnType computeType1Contribution(const PtBase<ScalarType>& a_ref_pt,
                                    const PtBase<ScalarType>& a_pt_0,
                                    const PtBase<ScalarType>& a_pt_1,
                                    bool* a_skip_first, bool a_is_arc,
                                    const NormalBase<ScalarType>& a_face_normal,
                                    const UnsignedIndex_t a_proj_dir) {
  using ReturnScalarType = typename ReturnType::value_type;
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    if (*a_skip_first) {
      *a_skip_first = false;
      return ReturnType::fromScalarConstant(ReturnScalarType(0));
    } else {
      return ReturnType::fromScalarConstant(ReturnScalarType(
          (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / ScalarType(6.0) *
          ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
           (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1]))));
    }
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    if (*a_skip_first) {
      *a_skip_first = false;
      return ReturnType::fromScalarConstant(ReturnScalarType(0));
    } else {
      /* Defining constants and types */
      const ScalarType ZERO = ScalarType(0);
      const ScalarType ONE = ScalarType(1);
      const ScalarType TWO = ScalarType(2);
      const ScalarType THREE = ScalarType(3);
      const ScalarType ONETWELVTH = ONE / ScalarType(12);

      /* Function */
      auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
      const ScalarType triangle_area =
          ((a_pt_0[0] - a_ref_pt[0]) * (a_pt_1[1] - a_ref_pt[1]) -
           (a_pt_1[0] - a_ref_pt[0]) * (a_pt_0[1] - a_ref_pt[1])) /
          TWO;
      moments.volume() = ReturnScalarType(
          triangle_area * (a_ref_pt[2] + a_pt_0[2] + a_pt_1[2]) / THREE);
      moments.centroid()[0] = ReturnScalarType(
          triangle_area *
          (a_pt_0[0] * (TWO * a_pt_0[2] + a_pt_1[2] + a_ref_pt[2]) +
           a_pt_1[0] * (a_pt_0[2] + TWO * a_pt_1[2] + a_ref_pt[2]) +
           a_ref_pt[0] * (a_pt_0[2] + a_pt_1[2] + TWO * a_ref_pt[2])) *
          ONETWELVTH);
      moments.centroid()[1] = ReturnScalarType(
          triangle_area *
          (a_pt_0[1] * (TWO * a_pt_0[2] + a_pt_1[2] + a_ref_pt[2]) +
           a_pt_1[1] * (a_pt_0[2] + TWO * a_pt_1[2] + a_ref_pt[2]) +
           a_ref_pt[1] * (a_pt_0[2] + a_pt_1[2] + TWO * a_ref_pt[2])) *
          ONETWELVTH);
      moments.centroid()[2] = ReturnScalarType(
          triangle_area *
          (a_pt_0[2] * a_pt_0[2] + a_pt_1[2] * a_pt_1[2] +
           a_ref_pt[2] * a_ref_pt[2] + a_pt_1[2] * a_ref_pt[2] +
           a_pt_0[2] * a_pt_1[2] + a_pt_0[2] * a_ref_pt[2]) *
          ONETWELVTH);
      return moments;
    }
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    if (!a_is_arc) {
      ParaboloidMomentLineIntegrator<ReturnType, ScalarType, 10> integrator(
          a_pt_0, a_pt_1, a_face_normal, a_proj_dir);
      return integrator.integrate();
    } else {
      return ReturnType::fromScalarConstant(ReturnScalarType(0));
    }
  } else {
    std::cout << "Type 1 for moments with order > 2 not yet implemented"
              << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

template <class ReturnType, class ScalarType>
ReturnType computeType2Contribution(
    const AlignedParaboloidBase<ScalarType>& a_aligned_paraboloid,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1) {
  using ReturnScalarType = typename ReturnType::value_type;
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    const ScalarType ONETWELVTH = ScalarType(1) / ScalarType(12);
    return ReturnType::fromScalarConstant(ReturnScalarType(
        ONETWELVTH * (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
        (-a_pt_0[2] - a_pt_1[2] +
         a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
         a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1])));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType ONETWELVTH = ONE / ScalarType(12);
    const ScalarType ONE60TH = ONE / ScalarType(60);
    const ScalarType ONE180TH = ONE / ScalarType(180);

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    moments.volume() = ReturnScalarType(
        (a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) * ONETWELVTH *
        (-a_pt_0[2] - a_pt_1[2] +
         a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] +
         a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1]));
    moments.centroid()[0] = ReturnScalarType(
        (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
        (TWO * a_aligned_paraboloid.b() * (a_pt_0[1] - a_pt_1[1]) *
             (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) +
         THREE * (a_pt_0[0] + a_pt_1[0]) * (a_pt_0[2] + a_pt_1[2])) *
        ONE60TH);
    moments.centroid()[1] = ReturnScalarType(
        (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
        (TWO * a_aligned_paraboloid.a() * (a_pt_0[0] - a_pt_1[0]) *
             (a_pt_1[1] * a_pt_0[0] - a_pt_0[1] * a_pt_1[0]) +
         THREE * (a_pt_0[1] + a_pt_1[1]) * (a_pt_0[2] + a_pt_1[2])) *
        ONE60TH);
    moments.centroid()[2] = ReturnScalarType(
        ((a_pt_0[0] * a_pt_1[1] - a_pt_1[0] * a_pt_0[1]) *
         (TWO * a_aligned_paraboloid.a() * a_aligned_paraboloid.b() *
              ((a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1]) *
               (a_pt_1[0] * a_pt_0[1] - a_pt_0[0] * a_pt_1[1])) +
          THREE * a_aligned_paraboloid.a() * a_pt_0[0] * a_pt_1[0] *
              (a_pt_0[2] + a_pt_1[2]) +
          THREE * a_aligned_paraboloid.b() * a_pt_0[1] * a_pt_1[1] *
              (a_pt_0[2] + a_pt_1[2]) -
          THREE * (a_pt_0[2] * a_pt_0[2] + a_pt_0[2] * a_pt_1[2] +
                   a_pt_1[2] * a_pt_1[2]))) *
        ONE180TH);
    return moments;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  } else {
    std::cout << "Type 2 for moments with order > 2 not yet implemented"
              << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsV3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 3> coeffs;
  coeffs.fill(ContainerType(ScalarType(0)));
  ContainerType x(ScalarType(1));
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      ContainerType add_to_coeff = ScalarType(v3Series[i][j]) * x;
      coeffs[j] += add_to_coeff;
      max_diff = maximum(max_diff, fabs(add_to_coeff));
    }
    if (max_diff < ScalarType(DBL_EPSILON)) {
      break;
    }
    x *= a_weight - ContainerType(ScalarType(1));
    i++;
  }
  return coeffs;
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 12> coeffsV3andC3SeriesOne(
    const ContainerType& a_weight) {
  std::array<ContainerType, 12> coeffs;
  coeffs.fill(ContainerType(ScalarType(0)));
  ContainerType x(ScalarType(1));
  UnsignedIndex_t i = 0;
  ScalarType max_diff;
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 3; ++j) {
      ContainerType add_to_coeff = ScalarType(v3Series[i][j]) * x;
      coeffs[j] += add_to_coeff;
      max_diff = maximum(max_diff, fabs(add_to_coeff));
    }
    if (max_diff < ScalarType(DBL_EPSILON)) {
      break;
    }
    x *= a_weight - ContainerType(ScalarType(1));
    i++;
  }
  i = 0;
  x = ContainerType(ScalarType(1));
  while (i <= 40) {
    max_diff = ScalarType(0);
    for (UnsignedIndex_t j = 0; j < 4; ++j) {
      ContainerType add_to_coeff = ScalarType(cx3Series[i][j]) * x;
      coeffs[3 + j] += add_to_coeff;
      max_diff = maximum(max_diff, fabs(add_to_coeff));
    }
    for (UnsignedIndex_t j = 0; j < 5; ++j) {
      ContainerType add_to_coeff = ScalarType(cz3Series[i][j]) * x;
      coeffs[7 + j] += add_to_coeff;
      max_diff = maximum(max_diff, fabs(add_to_coeff));
    }
    if (max_diff < ScalarType(DBL_EPSILON)) {
      break;
    }
    x *= a_weight - ContainerType(ScalarType(1));
    i++;
  }

  return coeffs;
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 3> coeffsV3Exact(
    const ContainerType& a_weight) {
  /* Defining constants and types */
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);
  const ScalarType FOUR = ScalarType(4);
  const ScalarType SIX = ScalarType(6);

  /* Function */
  const auto w2 = a_weight * a_weight;
  const auto w3 = w2 * a_weight;
  const auto w4 = w2 * w2;
  const auto w5 = w2 * w3;
  const auto w6 = w3 * w3;
  const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
                                       (a_weight + ContainerType(ONE)));
  const auto L3 = L * L * L;
  const auto S = (a_weight < ContainerType(ONE))
                     ? sqrt(ContainerType(ONE) - a_weight * a_weight)
                     : sqrt(a_weight * a_weight - ContainerType(ONE));
  const auto T = (a_weight < ContainerType(ONE))
                     ? atan((ContainerType(ONE) - a_weight) / S) / S
                     : atanh((a_weight - ContainerType(ONE)) / S) / S;
  return std::array<ContainerType, 3>(
      {(TWO * w6 - THREE * w4 + ScalarType(31) * w2 -
        (ScalarType(42) * w3 + ScalarType(18) * a_weight) * T) *
           L3 / ScalarType(12),
       (TWO * w6 - ScalarType(9) * w4 - ScalarType(8) * w2 +
        ScalarType(30) * w3 * T) *
           L3 / THREE,
       (ScalarType(11) * w4 + FOUR * w2 -
        (ScalarType(12) * w5 + ScalarType(18) * w3) * T) *
           L3 / SIX});
}

template <class ContainerType, class ScalarType>
inline std::array<ContainerType, 12> coeffsV3andC3Exact(
    const ContainerType& a_weight) {
  /* Defining constants and types */
  const ScalarType ONE = ScalarType(1);
  const ScalarType TWO = ScalarType(2);
  const ScalarType THREE = ScalarType(3);
  const ScalarType FOUR = ScalarType(4);
  const ScalarType FIVE = ScalarType(5);
  const ScalarType SIX = ScalarType(6);

  /* Function */
  const auto w2 = a_weight * a_weight;
  const auto w3 = w2 * a_weight;
  const auto w4 = w2 * w2;
  const auto w5 = w2 * w3;
  const auto w6 = w3 * w3;
  const auto w7 = w4 * w3;
  const auto w8 = w4 * w4;
  const auto w9 = w4 * w5;
  const auto w10 = w5 * w5;
  const auto L = ContainerType(ONE) / ((a_weight - ContainerType(ONE)) *
                                       (a_weight + ContainerType(ONE)));
  const auto L3 = L * L * L;
  const auto L4 = L3 * L;
  const auto L5 = L4 * L;
  const auto S = (a_weight < ContainerType(ONE))
                     ? sqrt(ContainerType(ONE) - a_weight * a_weight)
                     : sqrt(a_weight * a_weight - ContainerType(ONE));
  const auto T = (a_weight < ContainerType(ONE))
                     ? atan((ContainerType(ONE) - a_weight) / S) / S
                     : atanh((a_weight - ContainerType(ONE)) / S) / S;
  return std::array<ContainerType, 12>(
      {L3 *
           (TWO * w6 - THREE * w4 + ScalarType(31) * w2 -
            (ScalarType(42) * w3 + ScalarType(18) * a_weight) * T) /
           ScalarType(12),
       L3 *
           (TWO * w6 - ScalarType(9) * w4 - ScalarType(8) * w2 +
            ScalarType(30) * w3 * T) /
           THREE,
       L3 *
           (ScalarType(11) * w4 + FOUR * w2 -
            (ScalarType(12) * w5 + ScalarType(18) * w3) * T) /
           SIX,
       L4 * ((-T * a_weight) / ScalarType(32) +
             (ScalarType(93) * (w2)) / ScalarType(2240) -
             (ScalarType(163) * (w4)) / ScalarType(3360) +
             (FIVE * (w6)) / ScalarType(168) - (w8) / ScalarType(140)),
       L4 * ((w2) / ScalarType(70) + (-T * (w3)) / ScalarType(16) +
             (ScalarType(29) * (w4)) / ScalarType(1120) -
             (ScalarType(19) * (w6)) / ScalarType(1680) +
             (w8) / ScalarType(420)),
       -L4 *
           ((w2) / ScalarType(210) - (w4) / ScalarType(21) -
            (-T * (w5)) / ScalarType(8) -
            (ScalarType(13) * (w6)) / ScalarType(560) + (w8) / ScalarType(280)),
       L4 * ((w2) / ScalarType(35) - (ScalarType(16) * (w4)) / ScalarType(105) +
             (ScalarType(58) * (w6)) / ScalarType(105) - T * (w7) +
             (w8) / ScalarType(14)),
       L5 * ((-T * a_weight) / ScalarType(128) +
             (ScalarType(193) * (w2)) / ScalarType(16128) -
             (ScalarType(149) * (w4)) / ScalarType(8064) +
             (ScalarType(19) * (w6)) / ScalarType(1120) -
             (ScalarType(41) * (w8)) / ScalarType(5040) +
             (w10) / ScalarType(630)),
       L5 * ((FOUR * (w2)) / ScalarType(945) + (-T * (w3)) / ScalarType(48) +
             (ScalarType(65) * (w4)) / ScalarType(6048) -
             (w6) / ScalarType(144) +
             (ScalarType(11) * (w8)) / ScalarType(3780) -
             (w10) / ScalarType(1890)),
       -L5 * ((w2) / ScalarType(1890) -
              (ScalarType(13) * (w4)) / ScalarType(1890) -
              (-T * (w5)) / ScalarType(48) -
              (ScalarType(11) * (w6)) / ScalarType(2016) +
              (ScalarType(5) * (w8)) / ScalarType(3024) -
              (w10) / ScalarType(3780)),
       L5 * ((w2) / ScalarType(315) - (w4) / ScalarType(45) +
             (ScalarType(4) * (w6)) / ScalarType(35) +
             (-T * (w7)) / ScalarType(4) +
             (ScalarType(17) * (w8)) / ScalarType(504) -
             (w10) / ScalarType(252)),
       -L5 *
           ((w2) / ScalarType(63) - (ScalarType(29) * (w4)) / ScalarType(315) +
            (ScalarType(26) * (w6)) / ScalarType(105) -
            (ScalarType(194) * (w8)) / ScalarType(315) + T * (w9) -
            (w10) / ScalarType(18))});
}

template <class ReturnType, class ScalarType>
ReturnType computeType3Contribution(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const RationalBezierArcBase<ScalarType>& a_arc,
    const NormalBase<ScalarType>& a_face_normal) {
  using ReturnScalarType = typename ReturnType::value_type;
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);
    const ScalarType SIX = ScalarType(6);
    const ScalarType HALF = ONE / TWO;
    const ScalarType QUARTER = ONE / FOUR;

    const auto& pt_0 = a_arc.start_point();
    const auto& cp = a_arc.control_point();
    const auto& pt_1 = a_arc.end_point();
    const auto& weight = a_arc.weight();
    const ScalarType area_proj_triangle =
        HALF * (pt_0[0] * (pt_1[1] - cp[1]) + pt_1[0] * (cp[1] - pt_0[1]) +
                cp[0] * (pt_0[1] - pt_1[1]));
    assert(weight >= ZERO);
    std::array<ScalarType, 3> coeffs;
    if (weight < ScalarType(0.35))  // We use the exact expressions
      coeffs = coeffsV3Exact<ScalarType, ScalarType>(weight);
    else if (weight <
             ScalarType(1.7))  // We use the 40th order Taylor series (w -> 1)
      coeffs = coeffsV3SeriesOne<ScalarType, ScalarType>(weight);
    else if (weight <
             ScalarType(1.0e9))  // We use the series expansion (w -> infty)
      coeffs = coeffsV3Exact<ScalarType, ScalarType>(weight);
    else  // This is within EPSILON of the actual value
      coeffs = std::array<ScalarType, 3>({ONE / SIX, TWO / THREE, ZERO});

    return ReturnType::fromScalarConstant(ReturnScalarType(
        area_proj_triangle *
        (coeffs[0] *
             signedDistance<ScalarType>(HALF * (pt_0 + pt_1), a_paraboloid) +
         coeffs[1] * signedDistance<ScalarType>(
                         QUARTER * (pt_0 + pt_1) + HALF * cp, a_paraboloid) +
         coeffs[2] * signedDistance<ScalarType>(cp, a_paraboloid))));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);
    const ScalarType SIX = ScalarType(6);
    const ScalarType HALF = ONE / TWO;
    const ScalarType QUARTER = ONE / FOUR;

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    const auto& pt_0 = a_arc.start_point();
    const auto& cp = a_arc.control_point();
    const auto& pt_1 = a_arc.end_point();
    const auto weight = a_arc.weight();
    const auto A = a_paraboloid.a(), B = a_paraboloid.b();
    const auto X0 = pt_0[0], X1 = cp[0], X2 = pt_1[0];
    const auto Y0 = pt_0[1], Y1 = cp[1], Y2 = pt_1[1];
    const auto Z0 = pt_0[2], Z1p = cp[2], Z2 = pt_1[2];
    const ScalarType Z1P = -A * X1 * X1 - B * Y1 * Y1;
    const ScalarType AA = A * A, BB = B * B, AB = A * B;
    const ScalarType X00 = X0 * X0, X11 = X1 * X1, X22 = X2 * X2;
    const ScalarType X000 = X00 * X0, X111 = X11 * X1, X222 = X22 * X2;
    const ScalarType X0000 = X00 * X00, X1111 = X11 * X11, X2222 = X22 * X22;
    const ScalarType Y00 = Y0 * Y0, Y11 = Y1 * Y1, Y22 = Y2 * Y2;
    const ScalarType Y000 = Y00 * Y0, Y111 = Y11 * Y1, Y222 = Y22 * Y2;
    const ScalarType Y0000 = Y00 * Y00, Y1111 = Y11 * Y11, Y2222 = Y22 * Y22;
    const ScalarType Z00 = Z0 * Z0, Z22 = Z2 * Z2;
    const ScalarType Z1p1p = Z1p * Z1p, Z1P1P = Z1P * Z1P;
    const ScalarType X02 = X0 * X2, X12 = X1 * X2, X01 = X0 * X1;
    const ScalarType Y02 = Y0 * Y2, Y12 = Y1 * Y2, Y01 = Y0 * Y1;
    const ScalarType Z02 = Z0 * Z2, Z01p = Z0 * Z1p, Z01P = Z0 * Z1P;
    const ScalarType Z1p1P = Z1p * Z1P, Z1p2 = Z1p * Z2, Z1P2 = Z1P * Z2;
    const ScalarType X0Z0 = X0 * Z0, X0Z1p = X0 * Z1p, X0Z1P = X0 * Z1P,
                     X0Z2 = X0 * Z2;
    const ScalarType X1Z0 = X1 * Z0, X1Z1p = X1 * Z1p, X1Z1P = X1 * Z1P,
                     X1Z2 = X1 * Z2;
    const ScalarType X2Z0 = X2 * Z0, X2Z1p = X2 * Z1p, X2Z1P = X2 * Z1P,
                     X2Z2 = X2 * Z2;
    const ScalarType Y0Z0 = Y0 * Z0, Y0Z1p = Y0 * Z1p, Y0Z1P = Y0 * Z1P,
                     Y0Z2 = Y0 * Z2;
    const ScalarType Y1Z0 = Y1 * Z0, Y1Z1p = Y1 * Z1p, Y1Z1P = Y1 * Z1P,
                     Y1Z2 = Y1 * Z2;
    const ScalarType Y2Z0 = Y2 * Z0, Y2Z1p = Y2 * Z1p, Y2Z1P = Y2 * Z1P,
                     Y2Z2 = Y2 * Z2;
    const ScalarType area_proj_triangle =
        HALF * (X0 * (Y2 - Y1) + X2 * (Y1 - Y0) + X1 * (Y0 - Y2));
    assert(weight >= ZERO);
    // Compute coefficients (functions of the weight)
    std::array<ScalarType, 12> coeffs;
    if (weight < ScalarType(0.35))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
    } else if (weight <
               ScalarType(1.7))  // We use the 40th order Taylor series (w -> 1)
    {
      coeffs = coeffsV3andC3SeriesOne<ScalarType, ScalarType>(weight);
    } else if (weight < ScalarType(1.0e9))  // We use the exact expressions
    {
      coeffs = coeffsV3andC3Exact<ScalarType, ScalarType>(weight);
    }
    // else if (weight < 1.0e9)  // We use the series expansion (w -> infty)
    //   coeffs = coeffsV3andC3SeriesInfinity(weight);
    else  // This is within EPSILON of the actual value
    {
      coeffs = std::array<ScalarType, 12>(
          {ScalarType(ONE / SIX), ScalarType(TWO / THREE), ScalarType(ZERO),
           ScalarType(-ONE / ScalarType(140)),
           ScalarType(ONE / ScalarType(420)),
           ScalarType(-ONE / ScalarType(280)), ScalarType(ONE / ScalarType(14)),
           ScalarType(ONE / ScalarType(630)),
           ScalarType(-ONE / ScalarType(1890)),
           ScalarType(ONE / ScalarType(3780)),
           ScalarType(-ONE / ScalarType(252)),
           ScalarType(ONE / ScalarType(18))});
    }
    auto m0_basis = std::array<ScalarType, 3>(
        {signedDistance<ScalarType>(HALF * (pt_0 + pt_1), a_paraboloid),
         signedDistance<ScalarType>(QUARTER * (pt_0 + pt_1) + HALF * cp,
                                    a_paraboloid),
         signedDistance<ScalarType>(cp, a_paraboloid)});
    auto m1x_basis = std::array<ScalarType, 4>(
        {-SIX * (X0Z0 - X0Z2 - X2Z0 + X2Z2 - TWO * B * X2 * Y00 +
                 TWO * B * X0 * Y02 + TWO * B * X2 * Y02 - TWO * B * X0 * Y22),
         TWO * (FIVE * X0Z0 + ScalarType(10) * X0Z1p + SIX * X0Z1P +
                ScalarType(7) * X0Z2 + ScalarType(30) * A * X02 * X1 -
                ScalarType(11) * X1Z0 - FOUR * X1Z1p - ScalarType(11) * X1Z2 +
                ScalarType(7) * X2Z0 + ScalarType(10) * X2Z1p + SIX * X2Z1P +
                FIVE * X2Z2 - ScalarType(14) * B * X1 * Y00 +
                FOUR * B * X2 * Y00 + ScalarType(14) * B * X0 * Y01 -
                FOUR * B * X1 * Y01 + ScalarType(10) * B * X2 * Y01 -
                FOUR * B * X0 * Y02 + ScalarType(10) * B * X1 * Y02 -
                FOUR * B * X2 * Y02 + FOUR * B * X0 * Y11 +
                FOUR * B * X2 * Y11 + ScalarType(10) * B * X0 * Y12 -
                FOUR * B * X1 * Y12 + ScalarType(14) * B * X2 * Y12 +
                FOUR * B * X0 * Y22 - ScalarType(14) * B * X1 * Y22),
         TWO * (-FIVE * X0Z1p + ScalarType(18) * X0Z1P + X0Z2 +
                SIX * A * X02 * X1 - FIVE * X1Z0 - SIX * X1Z1p - SIX * X1Z1P -
                FIVE * X1Z2 + X2Z0 - FIVE * X2Z1p + ScalarType(18) * X2Z1P -
                ScalarType(12) * B * X1 * Y01 + TWO * B * X2 * Y01 +
                TWO * B * X1 * Y02 + ScalarType(12) * B * X0 * Y11 +
                ScalarType(12) * B * X2 * Y11 + TWO * B * X0 * Y12 -
                ScalarType(12) * B * X1 * Y12),
         TWO * (X1Z1p - X1Z1P)});
    auto m1y_basis = std::array<ScalarType, 4>(
        {SIX *
             (-Y0Z0 + Y0Z2 + TWO * A * (X22 * Y0 + X00 * Y2 - X02 * (Y0 + Y2)) +
              Y2Z0 - Y2Z2),
         TWO *
             (FIVE * Y0Z0 + ScalarType(10) * Y0Z1p + SIX * Y0Z1P +
              ScalarType(7) * Y0Z2 + ScalarType(30) * B * Y02 * Y1 -
              ScalarType(11) * Y1Z0 - FOUR * Y1Z1p - ScalarType(11) * Y1Z2 +
              TWO * A *
                  (-TWO * X02 * Y0 + TWO * X11 * Y0 + FIVE * X12 * Y0 +
                   TWO * X22 * Y0 - ScalarType(7) * X00 * Y1 + FIVE * X02 * Y1 -
                   TWO * X12 * Y1 - ScalarType(7) * X22 * Y1 + TWO * X00 * Y2 -
                   TWO * X02 * Y2 + TWO * X11 * Y2 + ScalarType(7) * X12 * Y2 +
                   X01 * (ScalarType(7) * Y0 - TWO * Y1 + FIVE * Y2)) +
              ScalarType(7) * Y2Z0 + ScalarType(10) * Y2Z1p + SIX * Y2Z1P +
              FIVE * Y2Z2),
         -TWO * (FIVE * Y0Z1p - ScalarType(18) * Y0Z1P - Y0Z2 -
                 SIX * B * Y02 * Y1 + FIVE * Y1Z0 + SIX * Y1Z1p + SIX * Y1Z1P +
                 FIVE * Y1Z2 -
                 TWO * A *
                     (X12 * Y0 - SIX * X01 * Y1 + X02 * Y1 - SIX * X12 * Y1 +
                      X01 * Y2 + SIX * X11 * (Y0 + Y2)) -
                 Y2Z0 + FIVE * Y2Z1p - ScalarType(18) * Y2Z1P),
         TWO * (Y1Z1p - Y1Z1P)});
    auto m1z_basis = std::array<ScalarType, 5>(
        {-(AA * (ScalarType(21) * X0000 + ScalarType(28) * X000 * X2 +
                 ScalarType(30) * X00 * X22 + ScalarType(28) * X0 * X222 +
                 ScalarType(21) * X2222)) -
             ScalarType(21) * BB * Y0000 - ScalarType(28) * BB * Y000 * Y2 -
             ScalarType(30) * BB * Y00 * Y22 -
             TWO * AB *
                 (X00 * (ScalarType(21) * Y00 + ScalarType(14) * Y02 +
                         FIVE * Y22) +
                  TWO * X02 *
                      (ScalarType(7) * Y00 + ScalarType(10) * Y02 +
                       ScalarType(7) * Y22) +
                  X22 * (ScalarType(5) * Y00 + ScalarType(14) * Y02 +
                         ScalarType(21) * Y22)) -
             ScalarType(28) * BB * Y0 * Y222 - ScalarType(21) * BB * Y2222 +
             ScalarType(40) * Z00 + ScalarType(48) * Z02 + ScalarType(40) * Z22,
         THREE * AA *
                 (ScalarType(21) * X000 * X1 - ScalarType(7) * X00 * X11 -
                  ScalarType(10) * X02 * X11 + ScalarType(35) * X00 * X12 -
                  ScalarType(7) * X000 * X2 - ScalarType(10) * X00 * X22 +
                  ScalarType(35) * X01 * X22 - ScalarType(7) * X11 * X22 -
                  ScalarType(7) * X0 * X222 + ScalarType(21) * X1 * X222) -
             AB * (ScalarType(7) * X11 * Y00 - ScalarType(35) * X12 * Y00 +
                   ScalarType(10) * X22 * Y00 - ScalarType(63) * X00 * Y01 +
                   ScalarType(20) * X12 * Y01 - ScalarType(35) * X22 * Y01 +
                   ScalarType(21) * X00 * Y02 + ScalarType(10) * X11 * Y02 -
                   ScalarType(70) * X12 * Y02 + ScalarType(21) * X22 * Y02 +
                   ScalarType(7) * X00 * Y11 + ScalarType(7) * X22 * Y11 -
                   ScalarType(35) * X00 * Y12 + ScalarType(28) * X12 * Y12 -
                   ScalarType(63) * X22 * Y12 +
                   X01 * (-ScalarType(63) * Y00 + ScalarType(28) * Y01 -
                          ScalarType(70) * Y02 + ScalarType(20) * Y12 -
                          ScalarType(35) * Y22) +
                   ScalarType(10) * X00 * Y22 + ScalarType(7) * X11 * Y22 -
                   ScalarType(63) * X12 * Y22 +
                   X02 * (ScalarType(21) * Y00 - ScalarType(70) * Y01 +
                          ScalarType(40) * Y02 + ScalarType(10) * Y11 -
                          ScalarType(70) * Y12 + ScalarType(21) * Y22)) -
             THREE *
                 (BB *
                      (ScalarType(7) * Y00 * Y11 + ScalarType(10) * Y02 * Y11 -
                       ScalarType(35) * Y00 * Y12 +
                       ScalarType(7) * Y000 * (-THREE * Y1 + Y2) +
                       ScalarType(10) * Y00 * Y22 - ScalarType(35) * Y01 * Y22 +
                       ScalarType(7) * Y11 * Y22 + ScalarType(7) * Y0 * Y222 -
                       ScalarType(21) * Y1 * Y222) +
                  TWO * (ScalarType(5) * Z00 + ScalarType(10) * Z01p +
                         FOUR * Z02 - TWO * Z1p1p + ScalarType(10) * Z1p2 +
                         FIVE * Z22)),
         -SIX * AA *
                 (ScalarType(46) * X02 * X11 - ScalarType(14) * X0 * X111 +
                  X1111 - ScalarType(14) * X111 * X2 -
                  ScalarType(14) * X01 * X22 + ScalarType(28) * X11 * X22 +
                  X00 * (ScalarType(28) * X11 - ScalarType(14) * X12 + X22)) -
             TWO * AB *
                 (X22 * Y00 + ScalarType(112) * X01 * Y01 -
                  ScalarType(28) * X02 * Y01 - ScalarType(14) * X22 * Y01 -
                  ScalarType(28) * X01 * Y02 + FOUR * X02 * Y02 +
                  ScalarType(28) * X00 * Y11 - ScalarType(42) * X01 * Y11 +
                  ScalarType(46) * X02 * Y11 + ScalarType(28) * X22 * Y11 -
                  TWO * X12 *
                      (ScalarType(7) * Y00 - ScalarType(46) * Y01 +
                       ScalarType(14) * Y02 + ScalarType(21) * Y11 -
                       ScalarType(56) * Y12) -
                  ScalarType(14) * X00 * Y12 + ScalarType(92) * X01 * Y12 -
                  ScalarType(28) * X02 * Y12 + X00 * Y22 -
                  ScalarType(14) * X01 * Y22 +
                  X11 * (ScalarType(28) * Y00 - ScalarType(42) * Y01 +
                         ScalarType(46) * Y02 + SIX * Y11 -
                         ScalarType(42) * Y12 + ScalarType(28) * Y22)) +
             THREE *
                 (-TWO * BB *
                      (ScalarType(46) * Y02 * Y11 - ScalarType(14) * Y0 * Y111 +
                       Y1111 - ScalarType(14) * Y111 * Y2 -
                       ScalarType(14) * Y01 * Y22 + ScalarType(28) * Y11 * Y22 +
                       Y00 * (ScalarType(28) * Y11 - ScalarType(14) * Y12 +
                              Y22)) +
                  FIVE * Z00 + ScalarType(40) * Z01p - TWO * Z02 +
                  ScalarType(8) * Z1p1p + ScalarType(40) * Z1p2 + FIVE * Z22),
         -TWO * AA *
                 (ScalarType(3) * X02 * X11 - ScalarType(7) * X0 * X111 +
                  THREE * X1111 - ScalarType(7) * X111 * X2) -
             SIX * BB * Y02 * Y11 + ScalarType(14) * BB * Y0 * Y111 -
             SIX * BB * Y1111 -
             TWO * AB *
                 (ScalarType(2) * X12 * Y01 - ScalarType(7) * X01 * Y11 +
                  X02 * Y11 - ScalarType(7) * X12 * Y11 +
                  X11 * (-ScalarType(7) * Y01 + Y02 + SIX * Y11 -
                         ScalarType(7) * Y12) +
                  TWO * X01 * Y12) +
             ScalarType(14) * BB * Y111 * Y2 - FIVE * Z01p + Z02 -
             ScalarType(7) * Z1p1p - FIVE * Z1p2,
         -(AA * X1111) - TWO * AB * X11 * Y11 - BB * Y1111 + Z1p1p});
    for (size_t i = 0; i < 3; ++i) {
      moments.volume() += ReturnScalarType(coeffs[i] * m0_basis[i]);
    }
    for (size_t i = 0; i < 4; ++i) {
      moments.centroid()[0] += ReturnScalarType(coeffs[3 + i] * m1x_basis[i]);
      moments.centroid()[1] += ReturnScalarType(coeffs[3 + i] * m1y_basis[i]);
    }
    for (size_t i = 0; i < 5; ++i) {
      moments.centroid()[2] += ReturnScalarType(coeffs[7 + i] * m1z_basis[i]);
    }
    moments.volume() *= ReturnScalarType(area_proj_triangle);
    moments.centroid()[0] *= ReturnScalarType(area_proj_triangle);
    moments.centroid()[1] *= ReturnScalarType(area_proj_triangle);
    moments.centroid()[2] *= ReturnScalarType(area_proj_triangle);
    return moments;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    ParaboloidMomentArcIntegrator<ReturnType, ScalarType, 10> integrator(
        a_paraboloid, a_arc, a_face_normal, 3);
    return integrator.integrate();
  } else {
    std::cout << "Type 3 for moments with order > 2 not yet implemented"
              << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

template <class ReturnType, class ScalarType>
ReturnType computeFaceOnlyContribution(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PlaneBase<ScalarType>& a_face_plane,
    const PtBase<ScalarType>& a_pt_ref) {
  using ReturnScalarType = typename ReturnType::value_type;
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    /* Defining constants and types */
    const ScalarType EPSILON = machine_epsilon<ScalarType>();
    const ScalarType ZERO = ScalarType(0);
    const ScalarType FOUR = ScalarType(4);

    /* Function */
    assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
    assert(fabs(a_face_plane.normal()[2]) > EPSILON);
    const ScalarType a =
        -a_face_plane.normal()[0] / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType b =
        -a_face_plane.normal()[1] / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType c =
        a_face_plane.distance() / safelyEpsilon(a_face_plane.normal()[2]);
    const ScalarType factor = FOUR * a_paraboloid.a() * a_paraboloid.b() * c -
                              a_paraboloid.a() * b * b -
                              a_paraboloid.b() * a * a;
    return ReturnType::fromScalarConstant(ReturnScalarType(copysign(
        machine_pi<ScalarType>() * factor * factor /
            (ScalarType(32) *
             pow(a_paraboloid.a() * a_paraboloid.b(), ScalarType(2.5))),
        -a_face_plane.normal()[2])));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    /* Defining constants and types */
    const ScalarType EPSILON = machine_epsilon<ScalarType>();
    const ScalarType ZERO = ScalarType(0);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);

    /* Function */
    assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
    assert(fabs(a_face_plane.normal()[2]) > EPSILON);
    const ScalarType a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
    const ScalarType b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
    const ScalarType c = a_face_plane.distance() / a_face_plane.normal()[2];
    const auto A = a_paraboloid.a(), B = a_paraboloid.b();
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    const ScalarType factor = (a * a * B + A * (b * b - FOUR * B * c)) *
                              (a * a * B + A * (b * b - FOUR * B * c)) *
                              machine_pi<ScalarType>();
    moments.volume() = ReturnScalarType(
        copysign(factor / (ScalarType(32) * pow(A * B, ScalarType(2.5))),
                 -a_face_plane.normal()[2]));
    moments.centroid()[0] = ReturnScalarType(
        a * B *
        copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]));
    moments.centroid()[1] = ReturnScalarType(
        b * A *
        copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]));
    moments.centroid()[2] = ReturnScalarType(
        (FIVE * A * (b * b) + FIVE * (a * a) * B - ScalarType(8) * A * B * c) *
        copysign(factor / (ScalarType(384) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]));
    return moments;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    /* Defining constants and types */
    const ScalarType EPSILON = machine_epsilon<ScalarType>();
    const ScalarType ZERO = ScalarType(0);
    const ScalarType TWO = ScalarType(2);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);
    const ScalarType SEVEN = ScalarType(7);

    /* Function */
    assert(a_paraboloid.a() * a_paraboloid.b() > ZERO);
    assert(fabs(a_face_plane.normal()[2]) > EPSILON);
    const ScalarType a = -a_face_plane.normal()[0] / a_face_plane.normal()[2];
    const ScalarType b = -a_face_plane.normal()[1] / a_face_plane.normal()[2];
    const ScalarType c = a_face_plane.distance() / a_face_plane.normal()[2];
    const auto A = a_paraboloid.a(), B = a_paraboloid.b();
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    const ScalarType factor = (a * a * B + A * (b * b - FOUR * B * c)) *
                              (a * a * B + A * (b * b - FOUR * B * c)) *
                              machine_pi<ScalarType>();
    moments[0] = ReturnScalarType(
        copysign(factor / (ScalarType(32) * pow(A * B, ScalarType(2.5))),
                 -a_face_plane.normal()[2]));
    moments[1] = ReturnScalarType(
        a * B *
        copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]));
    moments[2] = ReturnScalarType(
        b * A *
        copysign(factor / (ScalarType(64) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]));
    moments[3] = ReturnScalarType(
        (FIVE * A * (b * b) + FIVE * (a * a) * B - ScalarType(8) * A * B * c) *
        copysign(factor / (ScalarType(384) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]));
    moments[4] = ReturnScalarType(
        (A * B * b * b + SEVEN * B * B * a * a - FOUR * A * B * B * c) *
        copysign(factor / (ScalarType(768) * pow(A * B, ScalarType(4.5))),
                 a_face_plane.normal()[2]));
    moments[5] = ReturnScalarType(
        (a * b) *
        copysign(factor / (ScalarType(128) * pow(A * B, ScalarType(3.5))),
                 a_face_plane.normal()[2]));
    moments[6] = ReturnScalarType(
        (B * B * a * a * a + A * B * a * b * b - TWO * A * B * B * a * c) *
        copysign(factor / (ScalarType(128) * pow(A * B, ScalarType(4.5))),
                 a_face_plane.normal()[2]));
    moments[7] = ReturnScalarType(
        (A * B * a * a + SEVEN * A * A * b * b - FOUR * A * A * B * c) *
        copysign(factor / (ScalarType(768) * pow(A * B, ScalarType(4.5))),
                 a_face_plane.normal()[2]));
    moments[8] = ReturnScalarType(
        (A * A * b * b * b + A * B * a * a * b - TWO * A * A * B * b * c) *
        copysign(factor / (ScalarType(128) * pow(A * B, ScalarType(4.5))),
                 a_face_plane.normal()[2]));
    moments[9] = ReturnScalarType(
        (SEVEN * A * A * b * b * b * b + TWO * SEVEN * A * B * a * a * b * b +
         SEVEN * B * B * a * a * a * a -
         ScalarType(24) * A * A * B * b * b * c -
         ScalarType(24) * A * B * B * a * a * c +
         ScalarType(16) * A * A * B * B * c * c) *
        copysign(factor / (ScalarType(1024) * pow(A * B, ScalarType(4.5))),
                 a_face_plane.normal()[2]));
    return moments;
  } else {
    std::cout << "Type 4 for moments with order > 2 not yet implemented"
              << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

template <class ReturnType, class ScalarType>
ReturnType computeTriangleCorrection(
    const AlignedParaboloidBase<ScalarType>& a_paraboloid,
    const PtBase<ScalarType>& a_pt_0, const PtBase<ScalarType>& a_pt_1,
    const PtBase<ScalarType>& a_pt_2) {
  using ReturnScalarType = typename ReturnType::value_type;
  if constexpr (std::is_same_v<ReturnType, VolumeBase<ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(ReturnScalarType(
        (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
         -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
         a_pt_0[2] - ScalarType(2) * a_pt_1[2] - a_pt_2[2]) /
        ScalarType(12) *
        ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
         (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
         (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0])));
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  } else if constexpr (std::is_same_v<ReturnType,
                                      VolumeMomentsBase<ReturnScalarType>>) {
    /* Defining constants and types */
    const ScalarType ZERO = ScalarType(0);
    const ScalarType ONE = ScalarType(1);
    const ScalarType TWO = ScalarType(2);
    const ScalarType THREE = ScalarType(3);
    const ScalarType FOUR = ScalarType(4);
    const ScalarType FIVE = ScalarType(5);
    const ScalarType SIX = ScalarType(6);
    const ScalarType HALF = ONE / TWO;

    /* Function */
    auto moments = ReturnType::fromScalarConstant(ReturnScalarType(ZERO));
    const ScalarType A = a_paraboloid.a(), B = a_paraboloid.b();
    const ScalarType X0 = a_pt_0[0], X1 = a_pt_1[0], X2 = a_pt_2[0];
    const ScalarType Y0 = a_pt_0[1], Y1 = a_pt_1[1], Y2 = a_pt_2[1];
    const ScalarType Z0 = a_pt_0[2], Z1 = a_pt_1[2], Z2 = a_pt_2[2];
    const ScalarType triangle_area =
        HALF * ((a_pt_1[1] - a_pt_2[1]) * a_pt_0[0] +
                (a_pt_2[1] - a_pt_0[1]) * a_pt_1[0] +
                (a_pt_0[1] - a_pt_1[1]) * a_pt_2[0]);
    moments.volume() = ReturnScalarType(
        (-a_paraboloid.a() * (a_pt_0[0] + a_pt_1[0]) * (a_pt_1[0] + a_pt_2[0]) +
         -a_paraboloid.b() * (a_pt_0[1] + a_pt_1[1]) * (a_pt_1[1] + a_pt_2[1]) -
         a_pt_0[2] - TWO * a_pt_1[2] - a_pt_2[2]) *
        triangle_area / SIX);
    moments.centroid()[0] = ReturnScalarType(
        triangle_area *
        ((A * (-(X0 * X0 * X0) - X1 * X1 * X1 - X1 * X1 * X2 - X1 * (X2 * X2) -
               X2 * X2 * X2 - X0 * X0 * (X1 + X2) -
               X0 * (X1 * X1 + X1 * X2 + X2 * X2))) /
             ScalarType(10) +
         (B * (-(X1 * (Y0 * Y0 + TWO * Y0 * Y1 + THREE * (Y1 * Y1) + Y0 * Y2 +
                       TWO * Y1 * Y2 + Y2 * Y2)) -
               X2 * (Y0 * Y0 + Y0 * Y1 + Y1 * Y1 + TWO * Y0 * Y2 +
                     TWO * Y1 * Y2 + THREE * (Y2 * Y2)) -
               X0 * (THREE * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                     TWO * Y0 * (Y1 + Y2)))) /
             ScalarType(30) +
         (-(X0 * (TWO * Z0 + Z1 + Z2)) - X1 * (Z0 + TWO * Z1 + Z2) -
          X2 * (Z0 + Z1 + TWO * Z2)) /
             ScalarType(12)));
    moments.centroid()[1] = ReturnScalarType(
        -triangle_area *
        ((B * (Y0 * Y0 * Y0 + Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) +
               Y2 * Y2 * Y2 + Y0 * Y0 * (Y1 + Y2) +
               Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2))) /
             ScalarType(10) +
         (A *
          (X0 * X0 * (THREE * Y0 + Y1 + Y2) + X1 * X1 * (Y0 + THREE * Y1 + Y2) +
           X2 * X2 * (Y0 + Y1 + THREE * Y2) + X1 * X2 * (Y0 + TWO * (Y1 + Y2)) +
           X0 * (X1 * (TWO * Y0 + TWO * Y1 + Y2) +
                 X2 * (TWO * Y0 + Y1 + TWO * Y2)))) /
             ScalarType(30) +
         (Y0 * (TWO * Z0 + Z1 + Z2) + Y1 * (Z0 + TWO * Z1 + Z2) +
          Y2 * (Z0 + Z1 + TWO * Z2)) /
             ScalarType(12)));
    moments.centroid()[2] = ReturnScalarType(
        triangle_area *
        ((A * A *
          (X0 * X0 * X0 * X0 + X1 * X1 * X1 * X1 + X1 * X1 * X1 * X2 +
           X1 * X1 * (X2 * X2) + X1 * (X2 * X2 * X2) + X2 * X2 * X2 * X2 +
           X0 * X0 * X0 * (X1 + X2) + X0 * X0 * (X1 * X1 + X1 * X2 + X2 * X2) +
           X0 *
               (X1 * X1 * X1 + X1 * X1 * X2 + X1 * (X2 * X2) + X2 * X2 * X2))) /
             ScalarType(30) +
         (B * B *
          (Y0 * Y0 * Y0 * Y0 + Y1 * Y1 * Y1 * Y1 + Y1 * Y1 * Y1 * Y2 +
           Y1 * Y1 * (Y2 * Y2) + Y1 * (Y2 * Y2 * Y2) + Y2 * Y2 * Y2 * Y2 +
           Y0 * Y0 * Y0 * (Y1 + Y2) + Y0 * Y0 * (Y1 * Y1 + Y1 * Y2 + Y2 * Y2) +
           Y0 *
               (Y1 * Y1 * Y1 + Y1 * Y1 * Y2 + Y1 * (Y2 * Y2) + Y2 * Y2 * Y2))) /
             ScalarType(30) +
         (A * B *
          (X1 * X2 *
               (Y0 * Y0 + THREE * (Y1 * Y1) + FOUR * Y1 * Y2 +
                THREE * (Y2 * Y2) + TWO * Y0 * (Y1 + Y2)) +
           X0 * X0 *
               (SIX * (Y0 * Y0) + Y1 * Y1 + Y1 * Y2 + Y2 * Y2 +
                THREE * Y0 * (Y1 + Y2)) +
           X1 * X1 *
               (Y0 * Y0 + SIX * (Y1 * Y1) + THREE * Y1 * Y2 + Y2 * Y2 +
                Y0 * (THREE * Y1 + Y2)) +
           X2 * X2 *
               (Y0 * Y0 + Y1 * Y1 + THREE * Y1 * Y2 + SIX * (Y2 * Y2) +
                Y0 * (Y1 + THREE * Y2)) +
           X0 * (X1 * (THREE * (Y0 * Y0) + FOUR * Y0 * Y1 + THREE * (Y1 * Y1) +
                       TWO * Y0 * Y2 + TWO * Y1 * Y2 + Y2 * Y2) +
                 X2 * (THREE * (Y0 * Y0) + TWO * Y0 * Y1 + Y1 * Y1 +
                       FOUR * Y0 * Y2 + TWO * Y1 * Y2 + THREE * (Y2 * Y2))))) /
             ScalarType(90) +
         (-(Z0 * Z0) - Z1 * Z1 - Z1 * Z2 - Z2 * Z2 - Z0 * (Z1 + Z2)) /
             ScalarType(12)));
    return moments;
  } else if constexpr (std::is_same_v<
                           ReturnType,
                           GeneralMomentsBase<2, 3, ReturnScalarType>>) {
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  } else {
    std::cout << "Type 5 for moments with order > 2 not yet implemented"
              << std::endl;
    return ReturnType::fromScalarConstant(ReturnScalarType(0));
  }
}

}  // namespace IRL

#endif  // IRL_GENERIC_CUTTING_PARABOLOID_INTERSECTION_MOMENT_CONTRIBUTIONS_TPP_
