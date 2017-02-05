/**
 * @file Pdb Writer
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @private
 */
(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
    typeof define === 'function' && define.amd ? define(['exports'], factory) :
    (factory((global.Downloader = global.Downloader || {})));
}(this, (function (exports) { 'use strict';




var re = {
    not_string: /[^s]/,
    not_bool: /[^t]/,
    not_type: /[^T]/,
    not_primitive: /[^v]/,
    number: /[diefg]/,
    numeric_arg: /bcdiefguxX/,
    json: /[j]/,
    not_json: /[^j]/,
    text: /^[^\x25]+/,
    modulo: /^\x25{2}/,
    placeholder: /^\x25(?:([1-9]\d*)\$|\(([^\)]+)\))?(\+)?(0|'[^$])?(-)?(\d+)?(?:\.(\d+))?([b-gijostTuvxX])/,
    key: /^([a-z_][a-z_\d]*)/i,
    key_access: /^\.([a-z_][a-z_\d]*)/i,
    index_access: /^\[(\d+)\]/,
    sign: /^[\+\-]/
}

function sprintf() {
    var key = arguments[0], cache = sprintf.cache
    if (!(cache[key] && cache.hasOwnProperty(key))) {
        cache[key] = sprintf.parse(key)
    }
    return sprintf.format.call(null, cache[key], arguments)
}

sprintf.format = function(parse_tree, argv) {
    var cursor = 1, tree_length = parse_tree.length, node_type = '', arg, output = [], i, k, match, pad, pad_character, pad_length, is_positive = true, sign = ''
    for (i = 0; i < tree_length; i++) {
        node_type = get_type(parse_tree[i])
        if (node_type === 'string') {
            output[output.length] = parse_tree[i]
        }
        else if (node_type === 'array') {
            match = parse_tree[i] // convenience purposes only
            if (match[2]) { // keyword argument
                arg = argv[cursor]
                for (k = 0; k < match[2].length; k++) {
                    if (!arg.hasOwnProperty(match[2][k])) {
                        throw new Error(sprintf('[sprintf] property "%s" does not exist', match[2][k]))
                    }
                    arg = arg[match[2][k]]
                }
            }
            else if (match[1]) { // positional argument (explicit)
                arg = argv[match[1]]
            }
            else { // positional argument (implicit)
                arg = argv[cursor++]
            }

            if (re.not_type.test(match[8]) && re.not_primitive.test(match[8]) && get_type(arg) == 'function') {
                arg = arg()
            }

            if (re.numeric_arg.test(match[8]) && (get_type(arg) != 'number' && isNaN(arg))) {
                throw new TypeError(sprintf("[sprintf] expecting number but found %s", get_type(arg)))
            }

            if (re.number.test(match[8])) {
                is_positive = arg >= 0
            }

            switch (match[8]) {
                case 'b':
                    arg = parseInt(arg, 10).toString(2)
                break
                case 'c':
                    arg = String.fromCharCode(parseInt(arg, 10))
                break
                case 'd':
                case 'i':
                    arg = parseInt(arg, 10)
                break
                case 'j':
                    arg = JSON.stringify(arg, null, match[6] ? parseInt(match[6]) : 0)
                break
                case 'e':
                    arg = match[7] ? parseFloat(arg).toExponential(match[7]) : parseFloat(arg).toExponential()
                break
                case 'f':
                    arg = match[7] ? parseFloat(arg).toFixed(match[7]) : parseFloat(arg)
                break
                case 'g':
                    arg = match[7] ? parseFloat(arg).toPrecision(match[7]) : parseFloat(arg)
                break
                case 'o':
                    arg = arg.toString(8)
                break
                case 's':
                    arg = String(arg)
                    arg = (match[7] ? arg.substring(0, match[7]) : arg)
                break
                case 't':
                    arg = String(!!arg)
                    arg = (match[7] ? arg.substring(0, match[7]) : arg)
                break
                case 'T':
                    arg = get_type(arg)
                    arg = (match[7] ? arg.substring(0, match[7]) : arg)
                break
                case 'u':
                    arg = parseInt(arg, 10) >>> 0
                break
                case 'v':
                    arg = arg.valueOf()
                    arg = (match[7] ? arg.substring(0, match[7]) : arg)
                break
                case 'x':
                    arg = parseInt(arg, 10).toString(16)
                break
                case 'X':
                    arg = parseInt(arg, 10).toString(16).toUpperCase()
                break
            }
            if (re.json.test(match[8])) {
                output[output.length] = arg
            }
            else {
                if (re.number.test(match[8]) && (!is_positive || match[3])) {
                    sign = is_positive ? '+' : '-'
                    arg = arg.toString().replace(re.sign, '')
                }
                else {
                    sign = ''
                }
                pad_character = match[4] ? match[4] === '0' ? '0' : match[4].charAt(1) : ' '
                pad_length = match[6] - (sign + arg).length
                pad = match[6] ? (pad_length > 0 ? str_repeat(pad_character, pad_length) : '') : ''
                output[output.length] = match[5] ? sign + arg + pad : (pad_character === '0' ? sign + pad + arg : pad + sign + arg)
            }
        }
    }
    return output.join('')
}

sprintf.cache = {}

sprintf.parse = function(fmt) {
    var _fmt = fmt, match = [], parse_tree = [], arg_names = 0
    while (_fmt) {
        if ((match = re.text.exec(_fmt)) !== null) {
            parse_tree[parse_tree.length] = match[0]
        }
        else if ((match = re.modulo.exec(_fmt)) !== null) {
            parse_tree[parse_tree.length] = '%'
        }
        else if ((match = re.placeholder.exec(_fmt)) !== null) {
            if (match[2]) {
                arg_names |= 1
                var field_list = [], replacement_field = match[2], field_match = []
                if ((field_match = re.key.exec(replacement_field)) !== null) {
                    field_list[field_list.length] = field_match[1]
                    while ((replacement_field = replacement_field.substring(field_match[0].length)) !== '') {
                        if ((field_match = re.key_access.exec(replacement_field)) !== null) {
                            field_list[field_list.length] = field_match[1]
                        }
                        else if ((field_match = re.index_access.exec(replacement_field)) !== null) {
                            field_list[field_list.length] = field_match[1]
                        }
                        else {
                            throw new SyntaxError("[sprintf] failed to parse named argument key")
                        }
                    }
                }
                else {
                    throw new SyntaxError("[sprintf] failed to parse named argument key")
                }
                match[2] = field_list
            }
            else {
                arg_names |= 2
            }
            if (arg_names === 3) {
                throw new Error("[sprintf] mixing positional and named placeholders is not (yet) supported")
            }
            parse_tree[parse_tree.length] = match
        }
        else {
            throw new SyntaxError("[sprintf] unexpected placeholder")
        }
        _fmt = _fmt.substring(match[0].length)
    }
    return parse_tree
}

var vsprintf = function(fmt, argv, _argv) {
    _argv = (argv || []).slice(0)
    _argv.splice(0, 0, fmt)
    return sprintf.apply(null, _argv)
}

/**
 * helpers
 */
function get_type(variable) {
    if (typeof variable === 'number') {
        return 'number'
    }
    else if (typeof variable === 'string') {
        return 'string'
    }
    else {
        return Object.prototype.toString.call(variable).slice(8, -1).toLowerCase()
    }
}

var preformattedPadding = {
    '0': ['', '0', '00', '000', '0000', '00000', '000000', '0000000'],
    ' ': ['', ' ', '  ', '   ', '    ', '     ', '      ', '       '],
    '_': ['', '_', '__', '___', '____', '_____', '______', '_______'],
}
function str_repeat(input, multiplier) {
    if (multiplier >= 0 && multiplier <= 7 && preformattedPadding[input]) {
        return preformattedPadding[input][multiplier]
    }
    return Array(multiplier + 1).join(input)
}
function getBrowser(){

    if( typeof window === "undefined" ) return false;

    var ua = window.navigator.userAgent;

    if ( /Opera|OPR/.test( ua ) ) {

        return 'Opera';

    } else if ( /Chrome/i.test( ua ) ) {

        return 'Chrome';

    } else if ( /Firefox/i.test( ua ) ) {

        return 'Firefox';

    } else if ( /Mobile(\/.*)? Safari/i.test( ua ) ) {

        return 'Mobile Safari';

    } else if ( /MSIE/i.test( ua ) ) {

        return 'Internet Explorer';

    } else if ( /Safari/i.test( ua ) ) {

        return 'Safari';

    }

    return false;

}
function download( data, downloadName ){

    // using ideas from https://github.com/eligrey/FileSaver.js/blob/master/FileSaver.js

    if( !data ) return;

    downloadName = downloadName || "download";

    var isSafari = getBrowser() === "Safari";
    var isChromeIos = /CriOS\/[\d]+/.test( window.navigator.userAgent );

    var a = document.createElement( 'a' );

    function openUrl( url ){
        var opened = window.open( url, '_blank' );
        if( !opened ){
            window.location.href = url;
        }
    }

    function open( str ){
        openUrl( isChromeIos ? str : str.replace(/^data:[^;]*;/, 'data:attachment/file;') );
    }

    if( typeof navigator !== "undefined" && navigator.msSaveOrOpenBlob ){

        // native saveAs in IE 10+
        navigator.msSaveOrOpenBlob( data, downloadName );

    }else if( ( isSafari || isChromeIos ) && window.FileReader ){

        if( data instanceof Blob ){
            // no downloading of blob urls in Safari
            var reader = new FileReader();
            reader.onloadend = function() {
                open( reader.result );
            };
            reader.readAsDataURL( data );
        }else{
            open( data );
        }

    }else{

        if( data instanceof Blob ){
            data = URL.createObjectURL( data );
        }

        if( "download" in a ){
            // download link available
            a.style.display = "hidden";
            document.body.appendChild( a );
            a.href = data;
            a.download = downloadName;
            a.target = "_blank";
            a.click();
            document.body.removeChild( a );
        }else{
            openUrl( data );
        }

        if( data instanceof Blob ){
            URL.revokeObjectURL( data );
        }

    }

}
function LinkerPdbWriter( structure, linkerList, params ){


 






    var p = Object.assign( {}, params );

    var renumberSerial = p.renumberSerial !== undefined ? p.renumberSerial : true;
    var remarks = p.remarks || [];
    if( !Array.isArray( remarks ) ) remarks = [ remarks ];

    var records;

    function writeRecords(){

        records = [];

        writeTitle();
        writeRemarks();
        writeAtoms();

    }

    // http://www.wwpdb.org/documentation/file-format

    // Sample PDB line, the coords X,Y,Z are fields 5,6,7 on each line.
    // ATOM      1  N   ARG     1      29.292  13.212 -12.751  1.00 33.78      1BPT 108

    function DEF( x, y ){
        return x !== undefined ? x : y;
    }

    var atomFormat =
        "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s";

    var hetatmFormat =
        "HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s";

    function writeTitle(){

        // FIXME multiline if title line longer than 80 chars
        records.push( sprintf( "TITEL %-74s", structure.name ) );

    }

    function writeRemarks(){

        remarks.forEach( function( str ){
            records.push( sprintf( "REMARK %-73s", str ) );
        } );

        if( structure.trajectory ){
            records.push( sprintf(
                "REMARK %-73s",
                "Trajectory '" + structure.trajectory.name + "'"
            ) );
            records.push( sprintf(
                "REMARK %-73s",
                "Frame " + structure.trajectory.frame + ""
            ) );
        }

    }

    function writeAtoms(){

        var ia = 1;
        var im = 1;

        var linkerModels = [];
        var ntermResidues = [];
        var ntermChains = [];
        var ntermResnos = [];
        var ctermResidues = [];
        var ctermChains = [];
        var ctermResnos = [];

        for (var i = linkerList.length - 1; i >= 0; i--) {
            var linkerModel = linkerList[i];
            var linkerResidues = linkerModel.residueStore;
            var ntermResidue = linkerResidues.resno[ 0 ];
            ntermResidues.push(ntermResidue);
            var ctermResidue = linkerResidues.resno[ linkerResidues.count - 1 ];
            ctermResidues.push(ctermResidue);
            ntermChains.push(String.fromCharCode(linkerModel.chainStore.chainname[0]));//ntermResidue.atoms[ 0 ].chainname;
            ntermResnos.push(ntermResidue);
            ctermChains.push(String.fromCharCode(linkerModel.chainStore.chainname[0]));//ctermResidue.atoms[ 0 ].chainname;
            ctermResnos.push(ctermResidue);//.atoms[ 0 ].resno;
            linkerModels.push(linkerModel);
        }
        var currentCtermResno = 0;
        var inlinker = false;
        var modelc = 0;
        var linkerAdded = false;

        structure.eachModel( function( m ){

            records.push( sprintf( "MODEL %-74d", im++ ) );

            m.eachAtom( function( a ){
                for (var i = ntermResnos.length - 1; i >= 0; i--) {
                    if ( !linkerAdded && a.resno === ntermResnos[i]){
                        inlinker = true;
                        modelc = i;
                    }
                }
                if ( !linkerAdded && inlinker === true) {
                    linkerModels[modelc].eachAtom( function( b ){
                        currentCtermResno = ctermResnos[modelc];
                        var resno = b.resno;
                        if( resno >= ntermResnos[modelc] && resno <= ctermResnos[modelc] ){
                            var formatString = b.hetero ? hetatmFormat : atomFormat;
                            var serial = renumberSerial ? ia : b.serial;
                            var atomname = b.atomname;
                            if( atomname.length === 1 ) atomname = " " + atomname;

                            records.push( sprintf(
                                formatString,

                                serial,
                                atomname,
                                b.resname,
                                DEF( b.chainname, " " ),
                                b.resno,
                                b.x, b.y, b.z,
                                DEF( b.occurence, 1.0 ),
                                DEF( b.bfactor, 0.0 ),
                                DEF( b.segid, "" ),
                                DEF( b.element, "" )
                            ) );
                            ia += 1;
                        }
                    } );
                    linkerAdded = true;
                    inlinker = false;
                }else{
                    if (a.resno > currentCtermResno){ 
                        linkerAdded = false;
                        var formatString = a.hetero ? hetatmFormat : atomFormat;
                        var serial = renumberSerial ? ia : a.serial;

                        // Alignment of one-letter atom name such as C starts at column 14,
                        // while two-letter atom name such as FE starts at column 13.
                        var atomname = a.atomname;
                        if( atomname.length === 1 ) atomname = " " + atomname;

                        records.push( sprintf(
                            formatString,

                            serial,
                            atomname,
                            a.resname,
                            DEF( a.chainname, " " ),
                            a.resno,
                            a.x, a.y, a.z,
                            DEF( a.occurence, 1.0 ),
                            DEF( a.bfactor, 0.0 ),
                            DEF( a.segid, "" ),
                            DEF( a.element, "" )
                        ) );
                        ia += 1;
                    }
                }
            } );

            records.push( sprintf( "%-80s", "ENDMDL" ) );
            im += 1;

        } );

        records.push( sprintf( "%-80s", "END" ) );

    }

    function getString(){

        writeRecords();
        return records.join( "\n" );

    }

    function getBlob(){

        return new Blob(
            [ getString() ],
            { type: 'text/plain' }
        );

    }

    function _download( name, ext ){

        name = name || "structure";
        ext = ext || "pdb";

        var file = name + "." + ext;
        var blob = getBlob();

        download( blob, file );

    }

    // API

    this.getString = getString;
    this.getBlob = getBlob;
    this.download = _download;

}

function PdbWriter( structure, params ){

    var p = Object.assign( {}, params );

    var renumberSerial = p.renumberSerial !== undefined ? p.renumberSerial : true;
    var remarks = p.remarks || [];
    if( !Array.isArray( remarks ) ) remarks = [ remarks ];

    var records;

    function writeRecords(){

        records = [];

        writeTitle();
        writeRemarks();
        writeAtoms();

    }

    // http://www.wwpdb.org/documentation/file-format

    // Sample PDB line, the coords X,Y,Z are fields 5,6,7 on each line.
    // ATOM      1  N   ARG     1      29.292  13.212 -12.751  1.00 33.78      1BPT 108

    function DEF( x, y ){
        return x !== undefined ? x : y;
    }

    var atomFormat =
        "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s";

    var hetatmFormat =
        "HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s";

    function writeTitle(){

        // FIXME multiline if title line longer than 80 chars
        records.push( sprintf( "TITEL %-74s", structure.name ) );

    }

    function writeRemarks(){

        remarks.forEach( function( str ){
            records.push( sprintf( "REMARK %-73s", str ) );
        } );

        if( structure.trajectory ){
            records.push( sprintf(
                "REMARK %-73s",
                "Trajectory '" + structure.trajectory.name + "'"
            ) );
            records.push( sprintf(
                "REMARK %-73s",
                "Frame " + structure.trajectory.frame + ""
            ) );
        }

    }

    function writeAtoms(){

        var ia = 1;
        var im = 1;

        structure.eachModel( function( m ){

            records.push( sprintf( "MODEL %-74d", im++ ) );

            m.eachAtom( function( a ){

                var formatString = a.hetero ? hetatmFormat : atomFormat;
                var serial = renumberSerial ? ia : a.serial;

                // Alignment of one-letter atom name such as C starts at column 14,
                // while two-letter atom name such as FE starts at column 13.
                var atomname = a.atomname;
                if( atomname.length === 1 ) atomname = " " + atomname;

                records.push( sprintf(
                    formatString,

                    serial,
                    atomname,
                    a.resname,
                    DEF( a.chainname, " " ),
                    a.resno,
                    a.x, a.y, a.z,
                    DEF( a.occurence, 1.0 ),
                    DEF( a.bfactor, 0.0 ),
                    DEF( a.segid, "" ),
                    DEF( a.element, "" )
                ) );
                ia += 1;

            } );

            records.push( sprintf( "%-80s", "ENDMDL" ) );
            im += 1;

        } );

        records.push( sprintf( "%-80s", "END" ) );

    }

    function getString(){

        writeRecords();
        return records.join( "\n" );

    }

    function getBlob(){

        return new Blob(
            [ getString() ],
            { type: 'text/plain' }
        );

    }

    function _download( name, ext ){

        name = name || "structure";
        ext = ext || "pdb";

        var file = name + "." + ext;
        var blob = getBlob();

        download( blob, file );

    }

    // API

    this.getString = getString;
    this.getBlob = getBlob;
    this.download = _download;

}

exports.LinkerPdbWriter = LinkerPdbWriter;
exports.PdbWriter = PdbWriter;

Object.defineProperty(exports, '__esModule', { value: true });

})));
