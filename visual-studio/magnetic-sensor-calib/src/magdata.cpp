/**
  ******************************************************************************
  * @file    magdata.cpp
  * @author  @victorcroisfelt    
  * @version v2.0
  * @date    May 13, 2019.
  * @brief   Data to evaluate the calibration algorithms.
  ******************************************************************************
*/
/* Includes ------------------------------------------------------------------*/
#include <hdr\magdata.h>

/* Private types -------------------------------------------------------------*/

/* Private constants ---------------------------------------------------------*/

/* Private macro -------------------------------------------------------------*/

/* Private variables ---------------------------------------------------------*/

/* Private function prototypes -----------------------------------------------*/

/* Private functions ---------------------------------------------------------*/

/* Exported functions --------------------------------------------------------*/
/**
  ==============================================================================
  @def:				##### MagData Exported Functions #####
  ==============================================================================
  */
/**
  * @brief  Create the random data vectors.
  * @param  Return:'x','y','z', and the number of samples.
  * @retval None.
  */
void setupRawData(matrix **x,			
	matrix **y, 
	matrix **z,
	short sample_number) {
	/* Declaration of the three axis raw data (1 x N) */
	float axT[1][350] = { 
		206.468338012695312500f,
		200.419097900390625000f,
		187.336395263671875000f,
		164.818771362304687500f,
		119.653388977050781250f,
		64.209365844726562500f,
		-5.356473445892333984f,
		-95.066055297851562500f,
		-182.239929199218750000f,
		-266.618499755859375000f,
		-329.151123046875000000f,
		-353.367645263671875000f,
		-336.592407226562500000f,
		-291.204803466796875000f,
		-226.358810424804687500f,
		-154.166336059570312500f,
		-84.148437500000000000f,
		-28.019039154052734375f,
		11.174731254577636719f,
		27.662839889526367188f,
		17.850072860717773438f,
		-15.446266174316406250f,
		-78.105064392089843750f,
		-156.145736694335937500f,
		-243.006210327148437500f,
		-309.105560302734375000f,
		-338.504791259765625000f,
		-319.106781005859375000f,
		-254.931503295898437500f,
		-164.409179687500000000f,
		-64.298431396484375000f,
		23.766628265380859375f,
		99.865226745605468750f,
		145.225524902343750000f,
		170.974639892578125000f,
		164.514617919921875000f,
		131.960754394531250000f,
		71.871192932128906250f,
		-13.414504051208496094f,
		-111.468765258789062500f,
		-210.069961547851562500f,
		-293.034027099609375000f,
		-339.240325927734375000f,
		-334.965728759765625000f,
		-291.264190673828125000f,
		-212.533218383789062500f,
		-128.679779052734375000f,
		-54.581428527832031250f,
		-3.696682453155517578f,
		21.613218307495117188f,
		20.982780456542968750f,
		-4.490004062652587891f,
		-50.419647216796875000f,
		-114.377578735351562500f,
		-185.671508789062500000f,
		-261.949157714843750000f,
		-322.742340087890625000f,
		-357.235473632812500000f,
		-342.852264404296875000f,
		-281.481445312500000000f,
		-181.053573608398437500f,
		-60.764041900634765625f,
		52.902030944824218750f,
		138.425277709960937500f,
		188.854248046875000000f,
		203.976333618164062500f,
		186.383178710937500000f,
		144.037872314453125000f,
		85.147201538085937500f,
		14.193774223327636719f,
		-61.061923980712890625f,
		-132.596267700195312500f,
		-190.000152587890625000f,
		-221.282104492187500000f,
		-212.305953979492187500f,
		-172.554489135742187500f,
		-96.980224609375000000f,
		-11.854217529296875000f,
		61.593673706054687500f,
		108.360664367675781250f,
		115.760360717773437500f,
		92.849487304687500000f,
		44.993343353271484375f,
		-18.419670104980468750f,
		-88.662666320800781250f,
		-158.293304443359375000f,
		-224.921890258789062500f,
		-266.712921142578125000f,
		-285.991546630859375000f,
		-273.075286865234375000f,
		-220.692565917968750000f,
		-129.542190551757812500f,
		-11.549195289611816406f,
		109.132850646972656250f,
		207.817993164062500000f,
		264.832427978515625000f,
		276.836242675781250000f,
		247.324295043945312500f,
		189.685974121093750000f,
		119.237190246582031250f,
		48.980606079101562500f,
		-16.231510162353515625f,
		-66.972534179687500000f,
		-99.538856506347656250f,
		-105.817100524902343750f,
		-85.618385314941406250f,
		-36.017719268798828125f,
		38.535064697265625000f,
		123.873962402343750000f,
		200.829925537109375000f,
		247.058746337890625000f,
		251.054092407226562500f,
		204.209213256835937500f,
		130.037048339843750000f,
		38.252002716064453125f,
		-53.263141632080078125f,
		-130.983581542968750000f,
		-190.894149780273437500f,
		-228.551620483398437500f,
		-236.267898559570312500f,
		-216.925598144531250000f,
		-168.662902832031250000f,
		-91.548248291015625000f,
		7.425368785858154297f,
		115.830383300781250000f,
		214.215744018554687500f,
		283.114105224609375000f,
		299.482238769531250000f,
		269.569671630859375000f,
		203.742950439453125000f,
		120.723144531250000000f,
		43.787269592285156250f,
		-20.891290664672851563f,
		-61.512229919433593750f,
		-78.372383117675781250f,
		-68.328361511230468750f,
		-36.674438476562500000f,
		19.250068664550781250f,
		90.106315612792968750f,
		167.328735351562500000f,
		241.694366455078125000f,
		291.663177490234375000f,
		305.074188232421875000f,
		262.686645507812500000f,
		178.043960571289062500f,
		66.069725036621093750f,
		-45.012157440185546875f,
		-147.013916015625000000f,
		-213.065551757812500000f,
		-248.003387451171875000f,
		-247.622955322265625000f,
		-221.354248046875000000f,
		-172.742141723632812500f,
		-103.885955810546875000f,
		-25.808433532714843750f,
		55.480400085449218750f,
		130.593597412109375000f,
		180.070816040039062500f,
		198.906784057617187500f,
		175.837463378906250000f,
		115.259193420410156250f,
		33.920089721679687500f,
		-50.634346008300781250f,
		-112.921684265136718750f,
		-144.076995849609375000f,
		-142.823501586914062500f,
		-114.468406677246093750f,
		-64.144607543945312500f,
		-0.226012319326400757f,
		66.550788879394531250f,
		136.985809326171875000f,
		196.447982788085937500f,
		235.900299072265625000f,
		243.195800781250000000f,
		215.620437622070312500f,
		141.918426513671875000f,
		33.861888885498046875f,
		-90.060317993164062500f,
		-203.383453369140625000f,
		-281.253356933593750000f,
		-318.106750488281250000f,
		-310.896514892578125000f,
		-271.545867919921875000f,
		-211.923233032226562500f,
		-143.794174194335937500f,
		-73.431762695312500000f,
		-10.712699890136718750f,
		36.958438873291015625f,
		65.577705383300781250f,
		64.861717224121093750f,
		35.500705718994140625f,
		-24.678062438964843750f,
		-103.563331604003906250f,
		-190.945266723632812500f,
		-256.255920410156250000f,
		-283.301269531250000000f,
		-267.101074218750000000f,
		-213.963989257812500000f,
		-137.421295166015625000f,
		-50.460865020751953125f,
		31.690305709838867188f,
		100.929100036621093750f,
		152.071990966796875000f,
		179.621917724609375000f,
		180.975097656250000000f,
		151.328277587890625000f,
		88.060874938964843750f,
		-0.763267040252685547f,
		-108.022193908691406250f,
		-222.680175781250000000f,
		-312.285339355468750000f,
		-359.371368408203125000f,
		-357.162872314453125000f,
		-312.014007568359375000f,
		-238.870971679687500000f,
		-157.293106079101562500f,
		-83.874443054199218750f,
		-27.128906250000000000f,
		7.526398181915283203f,
		16.943258285522460938f,
		0.454189211130142212f,
		-37.908138275146484375f,
		-97.168792724609375000f,
		-173.364532470703125000f,
		-253.885177612304687500f,
		-323.157806396484375000f,
		-362.971069335937500000f,
		-352.227050781250000000f,
		-294.453948974609375000f,
		-198.274002075195312500f,
		-86.137229919433593750f,
		19.856603622436523438f,
		103.949272155761718750f,
		160.662353515625000000f,
		181.766128540039062500f,
		175.302902221679687500f,
		142.492279052734375000f,
		86.347045898437500000f,
		14.260489463806152344f,
		-72.460655212402343750f,
		-158.986740112304687500f,
		-234.959762573242187500f,
		-280.439727783203125000f,
		-287.744293212890625000f,
		-251.232696533203125000f,
		-177.618804931640625000f,
		-93.989822387695312500f,
		-15.287806510925292969f,
		37.480930328369140625f,
		64.636383056640625000f,
		57.550930023193359375f,
		26.830364227294921875f,
		-24.685668945312500000f,
		-90.655029296875000000f,
		-158.081481933593750000f,
		-226.699691772460937500f,
		-284.520904541015625000f,
		-319.818145751953125000f,
		-317.451049804687500000f,
		-274.284942626953125000f,
		-187.697875976562500000f,
		-69.238471984863281250f,
		52.862216949462890625f,
		154.310791015625000000f,
		220.185470581054687500f,
		243.195220947265625000f,
		227.658340454101562500f,
		184.330078125000000000f,
		124.183143615722656250f,
		54.767990112304687500f,
		-15.219936370849609375f,
		-77.352767944335937500f,
		-125.612762451171875000f,
		-151.929138183593750000f,
		-148.781860351562500000f,
		-111.577095031738281250f,
		-44.886787414550781250f,
		39.276336669921875000f,
		120.762031555175781250f,
		177.601287841796875000f,
		193.512115478515625000f,
		169.867660522460937500f,
		115.149398803710937500f,
		40.424343109130859375f,
		-39.391078948974609375f,
		-117.379959106445312500f,
		-178.805374145507812500f,
		-228.104690551757812500f,
		-252.754470825195312500f,
		-248.560440063476562500f,
		-211.550918579101562500f,
		-139.621658325195312500f,
		-37.734004974365234375f,
		78.822090148925781250f,
		189.744735717773437500f,
		267.911407470703125000f,
		298.673614501953125000f,
		291.376708984375000000f,
		233.638290405273437500f,
		160.652526855468750000f,
		82.731979370117187500f,
		12.739769935607910156f,
		-40.720260620117187500f,
		-73.091644287109375000f,
		-81.814720153808593750f,
		-63.812477111816406250f,
		-22.125402450561523438f,
		43.083129882812500000f,
		124.170501708984375000f,
		205.944137573242187500f,
		269.610443115234375000f,
		297.803741455078125000f,
		278.000427246093750000f,
		208.462036132812500000f,
		109.187675476074218750f,
		2.587900638580322266f,
		-95.179794311523437500f,
		-171.059478759765625000f,
		-218.807373046875000000f,
		-239.000427246093750000f,
		-227.768234252929687500f,
		-192.339599609375000000f,
		-131.116714477539062500f,
		-51.041095733642578125f,
		41.865631103515625000f,
		132.354370117187500000f,
		211.679779052734375000f,
		253.215866088867187500f,
		245.967742919921875000f,
		205.365859985351562500f,
		128.777603149414062500f,
		42.436008453369140625f,
		-32.866596221923828125f,
		-82.969345092773437500f,
		-106.828826904296875000f,
		-101.895240783691406250f,
		-70.470603942871093750f,
		-20.700984954833984375f,
		45.057781219482421875f,
		115.612457275390625000f,
		186.213668823242187500f,
		244.241302490234375000f,
		276.263000488281250000f,
		269.176940917968750000f,
		217.454071044921875000f,
		121.217697143554687500f,
		-2.218471765518188477f,
		-118.102317810058593750f,
		-211.901428222656250000f,
		-273.931091308593750000f
	};
		 
	float ayT[1][350] = {
		24.605636596679687500f,
		21.887414932250976563f,
		19.236852645874023438f,
		15.165770530700683594f,
		8.603707313537597656f,
		-0.516301095485687256f,
		-10.192257881164550781f,
		-30.515598297119140625f,
		-45.154800415039062500f,
		-63.234760284423828125f,
		-72.019462585449218750f,
		-71.369430541992187500f,
		-41.026145935058593750f,
		5.673428535461425781f,
		71.325256347656250000f,
		147.585464477539062500f,
		225.970962524414062500f,
		293.898101806640625000f,
		345.312591552734375000f,
		370.532470703125000000f,
		371.539703369140625000f,
		351.094421386718750000f,
		311.121124267578125000f,
		267.555023193359375000f,
		224.569107055664062500f,
		192.441452026367187500f,
		165.433685302734375000f,
		134.913955688476562500f,
		91.220405578613281250f,
		31.784528732299804688f,
		-40.772090911865234375f,
		-116.884284973144531250f,
		-186.321166992187500000f,
		-234.739501953125000000f,
		-257.668395996093750000f,
		-254.141616821289062500f,
		-219.598770141601562500f,
		-163.974258422851562500f,
		-92.222145080566406250f,
		-18.866212844848632813f,
		47.225372314453125000f,
		97.551681518554687500f,
		132.754501342773437500f,
		161.510070800781250000f,
		191.792083740234375000f,
		231.739669799804687500f,
		278.608703613281250000f,
		325.924377441406250000f,
		358.855560302734375000f,
		376.148590087890625000f,
		365.502044677734375000f,
		331.895111083984375000f,
		273.074859619140625000f,
		202.613769531250000000f,
		123.220909118652343750f,
		53.635112762451171875f,
		-4.035848617553710938f,
		-42.169227600097656250f,
		-67.681571960449218750f,
		-86.823287963867187500f,
		-109.315071105957031250f,
		-141.650588989257812500f,
		-179.452301025390625000f,
		-213.804382324218750000f,
		-232.879592895507812500f,
		-228.308853149414062500f,
		-196.366897583007812500f,
		-139.262756347656250000f,
		-61.384010314941406250f,
		27.392202377319335938f,
		116.655799865722656250f,
		199.262435913085937500f,
		260.341430664062500000f,
		302.993499755859375000f,
		321.630096435546875000f,
		329.395599365234375000f,
		333.535919189453125000f,
		339.774017333984375000f,
		346.153533935546875000f,
		346.742919921875000000f,
		332.596618652343750000f,
		296.208190917968750000f,
		238.676818847656250000f,
		158.888839721679687500f,
		72.039199829101562500f,
		-18.746446609497070313f,
		-96.677711486816406250f,
		-157.247680664062500000f,
		-193.713409423828125000f,
		-205.916641235351562500f,
		-199.228271484375000000f,
		-177.063323974609375000f,
		-158.029373168945312500f,
		-145.405853271484375000f,
		-140.341308593750000000f,
		-130.161819458007812500f,
		-110.330253601074218750f,
		-67.569953918457031250f,
		-4.217951297760009766f,
		75.833541870117187500f,
		158.141571044921875000f,
		237.789733886718750000f,
		305.498016357421875000f,
		351.085113525390625000f,
		370.430786132812500000f,
		365.782043457031250000f,
		342.074707031250000000f,
		308.617218017578125000f,
		273.565155029296875000f,
		242.145004272460937500f,
		220.772903442382812500f,
		193.454498291015625000f,
		156.354522705078125000f,
		100.334899902343750000f,
		27.156755447387695313f,
		-53.044921875000000000f,
		-131.294311523437500000f,
		-196.112976074218750000f,
		-237.680725097656250000f,
		-253.058914184570312500f,
		-239.684112548828125000f,
		-199.640258789062500000f,
		-145.249893188476562500f,
		-79.235275268554687500f,
		-21.362129211425781250f,
		26.482908248901367188f,
		62.314041137695312500f,
		89.779121398925781250f,
		120.920539855957031250f,
		165.326141357421875000f,
		219.460403442382812500f,
		276.910736083984375000f,
		328.842590332031250000f,
		364.441223144531250000f,
		377.002624511718750000f,
		361.585479736328125000f,
		323.130859375000000000f,
		264.178192138671875000f,
		195.624893188476562500f,
		125.104698181152343750f,
		63.304798126220703125f,
		17.662899017333984375f,
		-13.297504425048828125f,
		-36.936721801757812500f,
		-62.604232788085937500f,
		-100.143409729003906250f,
		-145.021423339843750000f,
		-191.385848999023437500f,
		-224.922409057617187500f,
		-240.763793945312500000f,
		-229.757293701171875000f,
		-190.474655151367187500f,
		-126.967071533203125000f,
		-48.152843475341796875f,
		39.175292968750000000f,
		122.120094299316406250f,
		193.351364135742187500f,
		245.831039428710937500f,
		276.965637207031250000f,
		293.829284667968750000f,
		308.537078857421875000f,
		324.974365234375000000f,
		342.868896484375000000f,
		357.013549804687500000f,
		359.515319824218750000f,
		340.837524414062500000f,
		297.743194580078125000f,
		232.777770996093750000f,
		150.948471069335937500f,
		61.707405090332031250f,
		-25.959659576416015625f,
		-97.817176818847656250f,
		-151.489227294921875000f,
		-179.787368774414062500f,
		-185.417083740234375000f,
		-176.877944946289062500f,
		-167.047210693359375000f,
		-160.528015136718750000f,
		-160.877822875976562500f,
		-162.277648925781250000f,
		-153.552108764648437500f,
		-124.999755859375000000f,
		-74.470573425292968750f,
		-4.665825843811035156f,
		78.898147583007812500f,
		163.951416015625000000f,
		245.874053955078125000f,
		307.429229736328125000f,
		348.508544921875000000f,
		362.903442382812500000f,
		362.623382568359375000f,
		339.809906005859375000f,
		316.149871826171875000f,
		294.261566162109375000f,
		275.685852050781250000f,
		256.796112060546875000f,
		225.663970947265625000f,
		177.332199096679687500f,
		109.765480041503906250f,
		26.962722778320312500f,
		-58.419281005859375000f,
		-140.808654785156250000f,
		-201.261337280273437500f,
		-239.339279174804687500f,
		-250.263946533203125000f,
		-233.349716186523437500f,
		-192.875152587890625000f,
		-139.170043945312500000f,
		-83.066085815429687500f,
		-34.176380157470703125f,
		1.034614205360412598f,
		28.044630050659179688f,
		53.590045928955078125f,
		96.996261596679687500f,
		149.583450317382812500f,
		212.653396606445312500f,
		276.080749511718750000f,
		329.088897705078125000f,
		366.256103515625000000f,
		375.875000000000000000f,
		360.556427001953125000f,
		322.222595214843750000f,
		269.358520507812500000f,
		205.614715576171875000f,
		146.887741088867187500f,
		96.992927551269531250f,
		61.369083404541015625f,
		32.084060668945312500f,
		2.808936595916748047f,
		-37.104381561279296875f,
		-88.438415527343750000f,
		-147.748184204101562500f,
		-200.939010620117187500f,
		-239.956069946289062500f,
		-252.281250000000000000f,
		-236.758819580078125000f,
		-192.813690185546875000f,
		-127.114807128906250000f,
		-47.268409729003906250f,
		37.693210601806640625f,
		116.884277343750000000f,
		179.983520507812500000f,
		223.586822509765625000f,
		253.300918579101562500f,
		269.552062988281250000f,
		293.470886230468750000f,
		318.194885253906250000f,
		342.039672851562500000f,
		362.067749023437500000f,
		363.981903076171875000f,
		346.566070556640625000f,
		300.417999267578125000f,
		233.582778930664062500f,
		154.104888916015625000f,
		69.158081054687500000f,
		-11.876796722412109375f,
		-76.663398742675781250f,
		-124.580329895019531250f,
		-148.744705200195312500f,
		-157.317016601562500000f,
		-155.173767089843750000f,
		-159.099975585937500000f,
		-168.590560913085937500f,
		-183.065139770507812500f,
		-191.025543212890625000f,
		-183.249740600585937500f,
		-153.923736572265625000f,
		-92.483428955078125000f,
		-16.284261703491210938f,
		70.163223266601562500f,
		158.487335205078125000f,
		239.098037719726562500f,
		302.752319335937500000f,
		342.204467773437500000f,
		359.034515380859375000f,
		353.626220703125000000f,
		340.668304443359375000f,
		323.762542724609375000f,
		309.896240234375000000f,
		297.981628417968750000f,
		279.386444091796875000f,
		243.754379272460937500f,
		190.635345458984375000f,
		117.569656372070312500f,
		33.186180114746093750f,
		-53.176620483398437500f,
		-131.110046386718750000f,
		-193.124191284179687500f,
		-227.817214965820312500f,
		-237.033477783203125000f,
		-220.892410278320312500f,
		-185.184097290039062500f,
		-141.800354003906250000f,
		-97.891052246093750000f,
		-65.694221496582031250f,
		-41.049873352050781250f,
		-17.961666107177734375f,
		14.795527458190917969f,
		60.727855682373046875f,
		125.505470275878906250f,
		196.303894042968750000f,
		266.846496582031250000f,
		324.971435546875000000f,
		363.159332275390625000f,
		376.125854492187500000f,
		366.597900390625000000f,
		330.386871337890625000f,
		277.147003173828125000f,
		220.943664550781250000f,
		168.297088623046875000f,
		126.577369689941406250f,
		96.353637695312500000f,
		66.862663269042968750f,
		30.388368606567382813f,
		-18.004581451416015625f,
		-78.347908020019531250f,
		-142.361145019531250000f,
		-199.552169799804687500f,
		-236.836074829101562500f,
		-249.802520751953125000f,
		-236.719818115234375000f,
		-194.210617065429687500f,
		-132.661819458007812500f,
		-55.943668365478515625f,
		21.378831863403320313f,
		94.152122497558593750f,
		149.283828735351562500f,
		187.734405517578125000f,
		215.457534790039062500f,
		240.792968750000000000f,
		269.377685546875000000f,
		305.303619384765625000f,
		339.752593994140625000f,
		364.407135009765625000f,
		370.320587158203125000f,
		351.916625976562500000f,
		310.283843994140625000f,
		246.882629394531250000f,
		165.737869262695312500f,
		83.881050109863281250f,
		6.762418270111083984f,
		-57.311283111572265625f,
		-102.833114624023437500f,
		-124.180358886718750000f,
		-137.654006958007812500f,
		-140.777877807617187500f,
		-155.640960693359375000f,
		-172.607330322265625000f,
		-194.607269287109375000f,
		-205.217910766601562500f
	};
	
		float azT[1][350] = { 
			294.748413085937500000f,
			299.740631103515625000f,
			313.838012695312500000f,
			337.462341308593750000f,
			362.207855224609375000f,
			382.427581787109375000f,
			396.735961914062500000f,
			390.441284179687500000f,
			349.550109863281250000f,
			273.934448242187500000f,
			164.707702636718750000f,
			36.515785217285156250f,
			-87.152137756347656250f,
			-185.936645507812500000f,
			-250.137313842773437500f,
			-271.111663818359375000f,
			-249.385330200195312500f,
			-193.220336914062500000f,
			-108.676872253417968750f,
			-5.321677207946777344f,
			97.838157653808593750f,
			187.525726318359375000f,
			247.265563964843750000f,
			267.012268066406250000f,
			231.022811889648437500f,
			145.431106567382812500f,
			22.689361572265625000f,
			-107.787132263183593750f,
			-219.992385864257812500f,
			-291.387023925781250000f,
			-312.182922363281250000f,
			-285.523559570312500000f,
			-216.069595336914062500f,
			-116.972190856933593750f,
			1.165558218955993652f,
			123.785659790039062500f,
			236.419906616210937500f,
			325.811035156250000000f,
			377.646240234375000000f,
			388.949127197265625000f,
			338.374938964843750000f,
			247.857025146484375000f,
			123.286331176757812500f,
			-7.349111080169677734f,
			-115.092941284179687500f,
			-177.463088989257812500f,
			-188.765838623046875000f,
			-150.130966186523437500f,
			-74.022773742675781250f,
			22.387605667114257813f,
			126.968025207519531250f,
			222.839065551757812500f,
			296.867797851562500000f,
			340.048156738281250000f,
			337.163146972656250000f,
			295.918029785156250000f,
			207.448196411132812500f,
			88.094749450683593750f,
			-47.630531311035156250f,
			-169.692962646484375000f,
			-254.292831420898437500f,
			-279.383728027343750000f,
			-242.423065185546875000f,
			-158.873718261718750000f,
			-43.894515991210937500f,
			80.953231811523437500f,
			200.148147583007812500f,
			299.443634033203125000f,
			366.743591308593750000f,
			396.979064941406250000f,
			387.685363769531250000f,
			336.006317138671875000f,
			250.280899047851562500f,
			139.701583862304687500f,
			28.660053253173828125f,
			-68.265495300292968750f,
			-122.895126342773437500f,
			-123.867576599121093750f,
			-71.608909606933593750f,
			19.254102706909179688f,
			127.264045715332031250f,
			235.153686523437500000f,
			319.206634521484375000f,
			374.044433593750000000f,
			390.911804199218750000f,
			367.392639160156250000f,
			305.898834228515625000f,
			211.474487304687500000f,
			92.931243896484375000f,
			-36.025470733642578125f,
			-152.769042968750000000f,
			-238.731002807617187500f,
			-269.516265869140625000f,
			-240.256317138671875000f,
			-152.400604248046875000f,
			-29.647188186645507813f,
			104.435859680175781250f,
			224.084518432617187500f,
			316.532196044921875000f,
			358.559814453125000000f,
			363.638031005859375000f,
			327.225982666015625000f,
			259.669342041015625000f,
			166.813064575195312500f,
			64.044052124023437500f,
			-37.942291259765625000f,
			-118.937736511230468750f,
			-160.126358032226562500f,
			-152.378128051757812500f,
			-89.111114501953125000f,
			16.374971389770507813f,
			141.622161865234375000f,
			259.280578613281250000f,
			348.727325439453125000f,
			392.962829589843750000f,
			389.732330322265625000f,
			341.565399169921875000f,
			257.466613769531250000f,
			150.702117919921875000f,
			28.466104507446289063f,
			-90.183982849121093750f,
			-196.575988769531250000f,
			-270.654754638671875000f,
			-301.701812744140625000f,
			-280.429351806640625000f,
			-204.696624755859375000f,
			-86.415092468261718750f,
			46.717666625976562500f,
			171.392333984375000000f,
			257.081298828125000000f,
			294.536315917968750000f,
			282.876342773437500000f,
			228.044418334960937500f,
			142.413024902343750000f,
			40.181869506835937500f,
			-60.757308959960937500f,
			-150.658401489257812500f,
			-216.379394531250000000f,
			-245.332824707031250000f,
			-230.621643066406250000f,
			-164.867385864257812500f,
			-58.516036987304687500f,
			72.969474792480468750f,
			205.897384643554687500f,
			309.753540039062500000f,
			360.980072021484375000f,
			353.772583007812500000f,
			293.110992431640625000f,
			194.444046020507812500f,
			75.327545166015625000f,
			-48.066497802734375000f,
			-160.916763305664062500f,
			-247.410629272460937500f,
			-304.710113525390625000f,
			-319.822448730468750000f,
			-291.707855224609375000f,
			-223.886108398437500000f,
			-122.931739807128906250f,
			-5.790967464447021484f,
			108.719230651855468750f,
			190.290359497070312500f,
			223.662704467773437500f,
			201.332656860351562500f,
			133.250610351562500000f,
			34.665565490722656250f,
			-72.809783935546875000f,
			-171.159225463867187500f,
			-248.461639404296875000f,
			-293.664947509765625000f,
			-304.035369873046875000f,
			-266.964752197265625000f,
			-193.709884643554687500f,
			-92.864341735839843750f,
			37.152477264404296875f,
			163.404052734375000000f,
			271.023376464843750000f,
			334.358673095703125000f,
			337.225036621093750000f,
			279.178558349609375000f,
			174.255981445312500000f,
			45.493129730224609375f,
			-82.102790832519531250f,
			-187.805908203125000000f,
			-262.103393554687500000f,
			-297.482757568359375000f,
			-287.078796386718750000f,
			-241.084259033203125000f,
			-163.350143432617187500f,
			-64.398651123046875000f,
			39.696254730224609375f,
			139.268951416015625000f,
			210.152435302734375000f,
			231.656646728515625000f,
			202.161956787109375000f,
			123.245010375976562500f,
			8.348496437072753906f,
			-113.133819580078125000f,
			-219.138977050781250000f,
			-290.521728515625000000f,
			-320.162170410156250000f,
			-303.822723388671875000f,
			-246.839355468750000000f,
			-158.860580444335937500f,
			-46.625144958496093750f,
			77.359031677246093750f,
			194.347015380859375000f,
			294.605499267578125000f,
			356.595916748046875000f,
			370.402404785156250000f,
			328.173126220703125000f,
			230.757186889648437500f,
			101.657554626464843750f,
			-37.148895263671875000f,
			-143.343215942382812500f,
			-213.353744506835937500f,
			-233.316711425781250000f,
			-208.148040771484375000f,
			-144.619964599609375000f,
			-51.729713439941406250f,
			50.297435760498046875f,
			152.386672973632812500f,
			238.073974609375000000f,
			293.419921875000000000f,
			310.545867919921875000f,
			278.853149414062500000f,
			196.571884155273437500f,
			77.233131408691406250f,
			-58.433792114257812500f,
			-181.553131103515625000f,
			-265.751464843750000000f,
			-296.389862060546875000f,
			-269.193847656250000000f,
			-198.387405395507812500f,
			-94.877044677734375000f,
			25.578477859497070313f,
			147.193435668945312500f,
			258.933197021484375000f,
			340.258148193359375000f,
			389.366943359375000000f,
			394.751770019531250000f,
			359.240142822265625000f,
			279.045135498046875000f,
			168.099624633789062500f,
			43.192115783691406250f,
			-63.990947723388671875f,
			-134.094390869140625000f,
			-150.487792968750000000f,
			-116.554550170898437500f,
			-40.070796966552734375f,
			61.159938812255859375f,
			166.073883056640625000f,
			261.492065429687500000f,
			331.042022705078125000f,
			369.855072021484375000f,
			370.145477294921875000f,
			328.289123535156250000f,
			245.095046997070312500f,
			132.475952148437500000f,
			2.132244586944580078f,
			-123.861518859863281250f,
			-219.889480590820312500f,
			-266.566192626953125000f,
			-249.025619506835937500f,
			-174.325195312500000000f,
			-61.476535797119140625f,
			70.602859497070312500f,
			190.862411499023437500f,
			289.956390380859375000f,
			357.396728515625000000f,
			385.768310546875000000f,
			373.869506835937500000f,
			324.121704101562500000f,
			242.697021484375000000f,
			141.354034423828125000f,
			33.686496734619140625f,
			-61.150131225585937500f,
			-124.006782531738281250f,
			-137.196960449218750000f,
			-92.947959899902343750f,
			-9.482174873352050781f,
			107.346710205078125000f,
			230.325225830078125000f,
			318.861724853515625000f,
			381.839752197265625000f,
			399.638397216796875000f,
			377.788543701171875000f,
			317.172515869140625000f,
			222.924530029296875000f,
			108.985015869140625000f,
			-14.539226531982421875f,
			-132.670028686523437500f,
			-225.174591064453125000f,
			-283.314208984375000000f,
			-276.700561523437500000f,
			-218.088119506835937500f,
			-107.544837951660156250f,
			26.793714523315429688f,
			157.418685913085937500f,
			258.372802734375000000f,
			319.736083984375000000f,
			332.416046142578125000f,
			301.698364257812500000f,
			234.486663818359375000f,
			143.778518676757812500f,
			37.588405609130859375f,
			-60.330741882324218750f,
			-146.423675537109375000f,
			-197.422393798828125000f,
			-208.260116577148437500f,
			-166.336166381835937500f,
			-77.231826782226562500f,
			47.850551605224609375f,
			180.862014770507812500f,
			294.174377441406250000f,
			364.258605957031250000f,
			381.971557617187500000f,
			347.636962890625000000f,
			270.867828369140625000f,
			163.109313964843750000f,
			43.355163574218750000f,
			-77.839080810546875000f,
			-182.400146484375000000f,
			-265.200256347656250000f,
			-310.285400390625000000f,
			-311.569213867187500000f,
			-265.249481201171875000f,
			-175.494644165039062500f,
			-55.662258148193359375f,
			70.861122131347656250f,
			182.066696166992187500f,
			240.105438232421875000f,
			250.447982788085937500f,
			208.184219360351562500f,
			131.170837402343750000f,
			31.441097259521484375f,
			-74.880790710449218750f,
			-168.293380737304687500f,
			-239.993637084960937500f,
			-279.092468261718750000f,
			-276.774749755859375000f,
			-233.287078857421875000f,
			-147.800125122070312500f,
			-31.498426437377929688f,
			102.869400024414062500f,
			225.162796020507812500f,
			314.960754394531250000f,
			348.957031250000000000f,
			321.374420166015625000f,
			239.581909179687500000f,
			124.311431884765625000f
		};
	
	    /* Attribution of 'x', 'y' and 'z' vectors (1 x N) */
	    (*x) = matrix2struct(axT[0],1,sample_number);
	    (*y) = matrix2struct(ayT[0],1,sample_number);
	    (*z) = matrix2struct(azT[0],1,sample_number);	
}
/******************************* END OF FILE **********************************/