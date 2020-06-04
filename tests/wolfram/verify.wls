file = OpenRead["output"];

f[irr_, prim_, m_, p_] :=
    irr == Boole[
        IrreduciblePolynomialQ[
            Dot[
                Power[
                    x,
                    Range[
                        0,
                        Length[p] - 1
                    ]
                ],
                p
            ],
            Modulus -> m
        ]
    ] &&
    prim == Boole[
        PrimitivePolynomialQ[
            Dot[
                Power[
                    x,
                    Range[
                        0,
                        Length[p] - 1
                    ]
                ],
                p
            ],
            Modulus -> m
        ]
    ];

While[
    Not[(irr = Read[file, Number]) === EndOfFile],
    prim = Read[file, Number];
    m = Read[file, Number];
    p = Read[StringToStream[#], Number] & /@
        StringSplit[
            StringTrim[
                ReadLine[file],
                ("{" | "}" | " ") ...
            ],
            ", "
        ];
    Print[f[irr, prim, m, p]]
];

Close[file];
