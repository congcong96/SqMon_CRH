function cmap = cschemes(map, nColors)
% cschemes  Color schemes for maps and graphics
%     cmap = cschemes(map, N) returns an RGB colormap specified 
%     by map with the number of colors specified by N. N is optional, and if
%     not supplied a default number of colors is used. Issuing colormap(cmap) 
%     in matlab then applies the desired colormap. Issuing cschemes without 
%     any input or output arguments displays and cycles through the possible 
%     colormaps using the matlab logo.
% 
%     There are two types of maps: diverging and sequential. Diverging is 
%     appropriate when values are centered around some meaningful, middle,
%     point. Sequential is best when values increase or decrease, and there 
%     is no logical division in the data.
% 
%     Diverging maps may be selected from the following options. Use the 
%     appropriate string as the value for map. The default number of colors
%     is 11 for diverging schems.
% 
%     1. 'puor' : purple/orange
%     2. 'brbg' : brown/blue-green
%     3. 'prgn' : purple/green
%     4. 'piyg' : purple/yellow-green
%     5. 'rdbu' : red/blue
%     6. 'rdgy' : red/gray
%     7. 'rdylbu' : red/yellow/blue
%     8. 'spectral' : red to green
%     9. 'rdylgn' : red/yellow/green
%     10. 'brbl' : brown/blue
%     11. 'bldo' : blue/dark orange
% 
%     For sequential maps, where the values go from some baseline upward
%     or downward, but not both, the following choices may be used. The 
%     default number for colors for sequential schemes is 9.
% 
%     1. 'ylgn' : yellow/green
%     2. 'ylgnbu': yellow-green-blue
%     3. 'gnbu': green/blue
%     4. 'bugn': blue/green
%     5. 'pubugn': purple/blue/green
%     6. 'pubu': purple/blue
%     7. 'bupu': blue/purple
%     8. 'rdpu': red/purple
%     9. 'purd': purple/red
%     10. 'orrd': orange/red
%     11. 'ylorrd': yellow/orange/red
%     12. 'ylorbr': yellow/orange/brown
%     13. 'purples': purples of different hues
%     14. 'blues'
%     15. 'greens'
%     16. 'oranges'
%     17. 'reds'
%     18. 'greys', 'grays'
% 
%     Colormaps are taken from the site http://ColorBrewer.org
% 
%     Craig Atencio
%     10/31/2011
% 
%     cmap = cschemes(map, N)

if ( nargin > 2 || nargout > 1 )
    error('Only 2 input args and 1 output arg allowed.');
end

divmaps = {'puor', 'brbg', 'prgn', 'piyg', 'rdbu', 'rdgy', 'rdylbu', ...
'spectral', 'rdylgn', 'brbl', 'bldo'};

seqmaps = {'ylgn', 'ylgnbu', 'gnbu', 'bugn', 'pubugn', 'pubu', ...
'bupu', 'rdpu', 'purd', 'orrd', 'ylorrd', 'ylorbr', 'purples', ...
'blues', 'greens', 'oranges', 'reds', 'greys', 'grays'};


if ( nargout == 0 && nargin == 0 ) % cycle through all the colormaps
    cyclemaps(divmaps, seqmaps);
end

if ( nargout == 1 && nargin == 0 ) % colormap not chosen
    error('You need to choose a colormap.');
end


if ( nargin == 1 || nargin == 2 ) % Get the colormap
    if ( ~ischar(map) )
        error('Map input arg must be a string.');
    end
    
    if ( sum( ismember(divmaps, map ) ) )
       cmap = diverging_maps(map);
    elseif ( sum( ismember(seqmaps, map ) ) )
       cmap = sequential_maps(map);
    else
       error('You chose an incorrect colormap.');
    end
end


if ( nargin == 2 ) % Calculate specified number of colors for colormap.
   cmap = resamplecmap(cmap, nColors);  
end


return;


% ------------------------------------------------------------------
% Function Definitions
% ------------------------------------------------------------------

function newcmap = resamplecmap(cmap, N)

   x = 1:size(cmap,1);
   newcmap = zeros(N,3);

   for i = 1:3 % go column by column
      y = cmap(:,i);
      xi = linspace(min(x), max(x), N);
      yi = interp1(x,y,xi,'linear');
      newcmap(:,i) = yi(:);
   end

   % Set any interpolation values < 0 or > 1 to be in [0,1]
   newcmap(newcmap>1) = 1;
   newcmap(newcmap<0) = 0;

return;


function cyclemaps(divmaps, seqmaps)

    figure;
    z = peaks;
    imagesc(z);
    colorbar;
    cmap = diverging_maps(divmaps{1});
    colormap(cmap);
    title('Colormaps')

    next_map = 1;
    i = 1;
    while ( next_map && i < length(divmaps) + 1 )
       cmap = diverging_maps(divmaps{i});
       colormap(cmap);
       title(sprintf('%s : %s', 'Diverging', divmaps{i}))
       i = i+1;
       next_map = input('Next (0 = stop, <enter> = continue) ? ');
       if isempty(next_map)
          next_map = 1;
       end
    end

    i = 1;
    while ( next_map && i < length(seqmaps) + 1 )
       cmap = sequential_maps(seqmaps{i});
       colormap(cmap);
       title(sprintf('%s : %s', 'Sequential', seqmaps{i}))
       i = i+1;
       next_map = input('Next (0 = stop, <enter> = continue) ? ');
       if isempty(next_map)
          next_map = 1;
       end
    end

    
return;



function cmap = diverging_maps(map)

if ( strcmp(map, 'puor') )

   x = [...
   127    59    8
   179    88    6
   224    130    20
   253    184    99
   254    224    182
   256    256    256
   216    218    235
   178    171    210
   128    115    172
   84    39    136
   45    0     75];

elseif ( strcmp(map, 'brbg') )

   x = [...
   84    48    5
   140    81    10
   191    129    45
   223    194    125
   246    232    195
   256    256    256
   199    234    229
   128    205    193
   53    151    143
   1     102    94
   0     60    48];

elseif ( strcmp(map, 'prgn') )

   x = [...
   64    0     75
   118    42    131
   153    112    171
   194    165    207
   231    212    232
   256    256    256
   217    240    211
   166    219    160
   90    174    97
   27    120    55
   0     68    27];

elseif ( strcmp(map, 'piyg') )

   x = [...
   142    1     82
   197    27    125
   222    119    174
   241    182    218
   253    224    239
   256    256    256
   230    245    208
   184    225    134
   127    188    65
   77    146    33
   39    100    25];

elseif ( strcmp(map, 'rdbu') )

   x = [...
   103    0     31
   178    24    43
   214   96    77
   244    165    130
   253    219    199
   256    256    256
   209    229    240
   146    197    222
   67    147    195
   33    102    172
   5     48    97];

elseif ( strcmp(map, 'rdgy') )

   x = [...
   103    0     31
   178    24    43
   214    96    77
   244    165    130
   253    219    199
   255    255    255
   224    224    224
   186    186    186
   135    135    135
   77    77    77
   26    26    26];

elseif ( strcmp(map, 'rdylbu') )

   x = [...
   165    0     38
   215    48    39
   244    109    67
   253    174    97
   254    224    144
   255    255    191
   224    243    248
   171    217    233
   116    173    209
   69    117    180
   49    54    149];

elseif ( strcmp(map, 'spectral') )

   x = [...
   158    1     66
   213    62    79
   244    109    67
   253    174    97
   254    224    139
   255    255    191
   230    245    152
   171    221    164
   102    194    165
   50    136    189
   94    79    162];

elseif ( strcmp(map, 'rdylgn') )

   x = [...
   165    0     38
   215    48    39
   244    109    67
   253    174    97
   254    224    139
   255    255    191
   217    239    139
   166    217    106
   102    189    99
   26    152    80
   0     104    55];


elseif ( strcmp(map, 'brbl') )

   x = [...
    51    25     0
   102    47     0
   153    96    53
   204   155   122
   216   175   151
   242   218   205
   256   256   256
   204   253   255
   153   248   255
   101   239   255
    50   227   255
     0   169   204
     0   122   153];


elseif ( strcmp(map, 'bldo') )

   x = [...
     0   102   102
     0   153   153
     0   204   204
     0   255   255
    51   255   255
   101   255   255
   153   255   255
   178   255   255
   203   255   255
   256   256   256
   229   255   255
   255   229   203
   255   202   153
   255   173   101
   255   142    51
   255   110     0
   204    85     0
   153    61     0
   102    39     0];

   x = flipud(x);

end

cmap = flipud(x)/256;

return;



function cmap = sequential_maps(map)


if ( strcmp(map, 'ylgn') )

   x = [...
   255	255	229
   247	252	185
   217	240	163
   173	221	142
   120	198	121
   65	171	93
   35	132	67
   0	104	55
   0	69	41];

elseif ( strcmp(map, 'ylgnbu') )

   x = [...
   255	255	217
   237	248	177
   199	233	180
   127	205	187
   65	182	196
   29	145	192
   34	94	168
   37	52	148
   8	29	88];

elseif ( strcmp(map, 'gnbu') )

   x = [...
   247	252	240
   224	243	219
   204	235	197
   168	221	181
   123	204	196
   78	179	211
   43	140	190
   8	104	172
   8	64	129];

elseif ( strcmp(map, 'bugn') )

   x = [...
   247	252	253
   229	245	249
   204	236	230
   153	216	201
   102	194	164
   65	174	118
   35	139	69
   0	109	44
   0	68	27];

elseif ( strcmp(map, 'pubugn') )

   x = [...
   255	247	251
   236	226	240
   208	209	230
   166	189	219
   103	169	207
   54	144	192
   2	129	138
   1	108	89
   1	70	54];

elseif ( strcmp(map, 'pubu') )

   x = [...
   255	247	251
   236	231	242
   208	209	230
   166	189	219
   116	169	207
   54	144	192
   5	112	176
   4	90	141
   2	56	88];

elseif ( strcmp(map, 'bupu') )

   x = [...
   247	252	253
   224	236	244
   191	211	230
   158	188	218
   140	150	198
   140	107	177
   136	65	157
   129	15	124
   77	0	75];

elseif ( strcmp(map, 'rdpu') )

   x = [...
   255	247	243
   253	224	221
   252	197	192
   250	159	181
   247	104	161
   221	52	151
   174	1	126
   122	1	119
   73	0	106];

elseif ( strcmp(map, 'purd') )

   x = [...
   247	244	249
   231	225	239
   212	185	218
   201	148	199
   223	101	176
   231	41	138
   206	18	86
   152	0	67
   103	0	31];

elseif ( strcmp(map, 'orrd') )

   x = [...
   255	247	236
   254	232	200
   253	212	158
   253	187	132
   252	141	89
   239	101	72
   215	48	31
   179	0	0
   127	0	0];

elseif ( strcmp(map, 'ylorrd') )

   x = [...
   255	255	204
   255	237	160
   254	217	118
   254	178	76
   253	141	60
   252	78	42
   227	26	28
   189	0	38
   128	0	38];

elseif ( strcmp(map, 'ylorbr') )

   x = [...
   255	255	229
   255	247	188
   254	227	145
   254	196	79
   254	153	41
   236	112	20
   204	76	2
   153	52	4
   102	37	6];

elseif ( strcmp(map, 'purples') )

   x = [...
   252	251	253
   239	237	245
   218	218	235
   188	189	220
   158	154	200
   128	125	186
   106	81	163
   84	39	143
   63	0	125];

elseif ( strcmp(map, 'blues') )

   x = [...
   247	251	255
   222	235	247
   198	219	239
   158	202	225
   107	174	214
   66	146	198
   33	113	181
   8	81	156
   8	48	107];

elseif ( strcmp(map, 'greens') )

   x = [...
   247	252	245
   229	245	224
   199	233	192
   161	217	155
   116	196	118
   65	171	93
   35	139	69
   0	109	44
   0	68	27];

elseif ( strcmp(map, 'oranges') )

   x = [...
   255	245	235
   254	230	206
   253	208	162
   253	174	107
   253	141	60
   241	105	19
   217	72	1
   166	54	3
   127	39	4];

elseif ( strcmp(map, 'reds') )

   x = [...
   255	245	240
   254	224	210
   252	187	161
   252	146	114
   251	106	74
   239	59	44
   203	24	29
   165	15	21
   103	0	13];

elseif ( strcmp(map, 'greys') ||  strcmp(map, 'grays') )

   x = [...
   255	255	255
   240	240	240
   217	217	217
   189	189	189
   150	150	150
   115	115	115
   82	82	82
   37	37	37
   0	0	0];

end

cmap = flipud(x)/256;


return;
















