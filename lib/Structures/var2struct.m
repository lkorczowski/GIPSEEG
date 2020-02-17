function varargout = var2struct(varargin)
% V2STRUCT Pack/Unpack Variables to/from a scalar structure.
%
%% Description:
%    V2STRUCT has dual functionality in packing & unpacking variables into structures and
%    vice versa, according to the syntax and inputs.
%
%    Function features:
%      * Pack variables to new structure with enhanced field naming
%      * Pack and update variables in existing structure
%      * Unpack variables from structure with enhanced variable naming
%      * Unpack only specific fields in a structure to variables
%      * Unpack without over writing existing variables in work space
%
%    In addition to the obvious usage, this function could by highly useful for example in
%    working with a function with multiple inputs. Packing variables before the call to
%    the function, and unpacking it in the beginning of the function will make the function
%    call shorter, more readable, and you would not have to worry about arguments order any
%    more. Moreover you could leave the function as it is and you could pass same inputs to
%    multiple functions, each of which will use its designated arguments placed in the
%    structure.
%
%% Syntax:
%    Pack:
%      S = V2STRUCT(x,y,z,...)
%      S = V2STRUCT(fieldNames)
%      S = V2STRUCT(A,B,C,...,fieldNames)
%      S = V2STRUCT(x,...,nameOfStruct2Update,fieldNames)
%      V2STRUCT(x,y,z,...)
%      V2STRUCT(fieldNames)
%      V2STRUCT(A,B,C,...,fieldNames)
%      V2STRUCT(x,...,nameOfStruct2Update,fieldNames)
%
%    Unpack:
%      V2STRUCT(S)
%      [a,b,c,...] = V2STRUCT(S)
%      V2STRUCT(S,fieldNames)
%      [a,b,c,...] = V2STRUCT(S,fieldNames)
%
%% Inputs & Outputs:
%    Pack - inputs:
%      (reminder: S = V2STRUCT(x,y,z,...,nameOfStruct2Update,fieldNames) )
%       x,y,z,... - any variable to pack. can be replaced by fieldNames below.
%       nameOfStruct2Update - optional, name of structure to update if desired.
%       fieldNames - optional, cell array of strings, which must include a cell with the string 'fieldNames' and must be the last input.
%    Pack - outputs: 
%       S - the packed structure
%       if there is no output argument then a structure named Sv2struct would be created in the caller workspace
%
%    Unpack - inputs:
%      (reminder: [a,b,c,...] = V2STRUCT(S,fieldNames) )
%       S - structure to be unpacked.
%       fieldNames - optional, cell array of strings, which must include a cell with the string 'fieldNames' and must be the last input.
%    Unpack - outputs:          
%       a,b,c,... - variables upacked from the structure
%       if there are no output arguments then variables would be created in the caller workspace with naming according to inputs.
%
%% Examples:
%  % Note: see Usage example section below for convenient presentation of these examples.
%
%    % NOTE: whenever using filedNames cell array please note the following -
%    %    1. fieldNames cell array must include a cell with the string 'fieldNames'
%    %    2. fieldNames cell array input must be the last input.
%
%  % Pack:
%      x = zeros(3); x2 = ones(3); y = 'Testing123'; z = cell(2,3);
%      fieldNames1 = {'fieldNames','x','y','z'};
%      fieldNames2 = {'fieldNames','a','b','c'};
%      fieldNames3 = {'fieldNames','x'};
%      nameOfStruct2Update = 'S';
%
%       % The four examples below return structure S with same values however the
%       %  structure'sfield names are defined differently in every syntax. 
%      % Example 1.
%      S = v2struct(x,y,z) % structure field names defined by variables names.
%      % Example 2.
%      S = v2struct(fieldNames1) % structure field names defined according to the cell array
%       % fieldNames. 
%       % NOTE: variables with the names in fieldNames1 must exist in the caller workspace.
%      % Example 3.
%      S = v2struct(zeros(3),'Testing123',cell(2,3),fieldNames1) % same as #1. but
%       % arguments are passed explicitly
%      % Example 4.
%      S = v2struct(x,y,z,fieldNames2) % field names defined by content of fieldNames2 while
%       % the values are set according to the passed arguments. In this case the structure
%       % S returned would be: S.a=x, S.b=y, S.c=z
%
%      % Example 5.
%      S.oldField = 'field to be saved for future use'
%      S = v2struct(x2,nameOfStruct2Update,fieldNames3) % update structure S. The fields
%       % that would be updated are according to content of fieldNames3. Note that you must
%       % pass a variable with the name 'nameOfStruct2Update' placed before 'fieldNames3'.
%       % This variable should contain the name of the structure you want to update as a
%       % string. Also note that if you set an output structure name which is different
%       % than the one stated in nameOfStruct2Update a new structure would be created and
%       % the structure that was meant to be updated would not get updated.
%
%       % The following examples return the same results as the examples above but the
%       % structure would be returned with the default name 'Sv2struct'. This might lead
%       % to overriding of arguments so this is unrecommended.
%      % Example 6.
%      v2struct(x,y,z)
%      % Example 7.
%      v2struct(fieldNames1)
%      % Example 8.
%      v2struct(zeros(3),'Testing123',cell(2,3),fieldNames1)
%      % Example 9.
%      v2struct(x,y,z,fieldNames2)
%      % Example 10.
%      S.oldField = 'field to be saved for future use'
%      v2struct(x2,nameOfStruct2Update,fieldNames3)
%
%  % Unpack:
%      S.x = zeros(3); S.y = 'Testing123'; S.z = cell(2,3);
%      fieldNames3 = {'fieldNames','x','z'};
%
%       % This example creates or overwrites variables x, y, z in the caller with the
%       % contents of the corresponding named fields.
%      % Example 1.
%      v2struct(S)
%
%       % This example assigns the contents of the fields of the scalar structure
%       % S to the variables a,b,c rather than overwriting variables in the caller. If
%       % there are fewer output variables than there are fields in S, the remaining fields
%       % are not extracted.
%      % Example 2.
%      [a,b,c] = v2struct(S)
%
%       % This example creates or overwrites variables x and z in the caller with the
%       % contents of the corresponding named fields.
%      % Example 3.
%      v2struct(S,fieldNames3)
%
%       % This example assigns the contents of the fields 'x' and 'z' defined by
%       % fieldNames3 of the scalar structure S to the variables a and b rather than
%       % overwriting variables in the caller. If there are fewer output variables than
%       % there are fields in S, the remaining fields are not extracted.
%      % Example 4.
%      [a,b] = v2struct(S,fieldNames3)
%
%       % This example unpacks variable 'z' only without overwriting variable 'x'. 
%       % NOTE: the addition of the field named 'avoidOverWrite' to the structure to be
%       % unpacked. This is mandatory in order to make this functionality work. The
%       % contents of this field can be anything, it does not matter. 
%      S.avoidOverWrite = 'foo(contents does not matter)';
%      x = 'do not overwrite me';
%      v2struct(S)
%
%% Usage example (includes sub-functions):
%    1. run attached v2structDemo1.m file for on screen presentation of examples.
%    2. run attached v2structDemo2.m file and read comments in file for a suggestion of
%       how to use v2struct in managing input to other functions with improved usability.
%
%% Revision history:
%    2011-05-19, Adi N., Creation
%    2011-05-29, Adi N., update structure added, some documentation and demo function changes
%    2011-06-02, Adi N., fixed updating structure functionality
%    2011-06-05, Adi N., Added functionality: avoid overwritring existing variables, added
%                        unpacking examples to demo1 .m file.
%    2011-06-30, Adi N., fieldNames usage corrected, now must include a specific string to
%                        be triggered. Documentation enhanced. Code tweaked.
%    2011-07-14, Adi N., fixed bug in packing with variables only
%
%    Inspired by the function: mmv2truct - D.C. Hanselman, University of Maine, Orono, ME
%    04469 4/28/99, 9/29/99, renamed 10/19/99 Mastering MATLAB 5, Prentice Hall,
%    ISBN 0-13-858366-8
%%

% parse input for field names
gotCellArrayOfStrings = iscellstr(varargin{end});

toUnpackRegular = nargin == 1 && isstruct(varargin{1});
if toUnpackRegular
   fieldNames = fieldnames(varargin{1})';
   nFields = length(fieldNames);
end

gotFieldNames = gotCellArrayOfStrings & any(strcmpi(varargin{end},'fieldNames'));
if gotFieldNames
   fieldNamesRaw = varargin{end};
   % indices of cells with actual field names, leaving out the index to 'fieldNames' cell
   indFieldNames = ~strcmpi(fieldNamesRaw,'fieldNames');
   fieldNames = fieldNamesRaw(indFieldNames);
   nFields = length(fieldNames);
end
toUnpackFieldNames = nargin == 2 && isstruct(varargin{1}) && gotFieldNames;


if nargin == 0
   error('v2struct: Input arguments required.')
   
elseif toUnpackRegular || toUnpackFieldNames % Unpack
   
   struct = varargin{1};
   assert(isequal(length(struct),1) , 'v2struct: Single input must be a scalar structure.');
   callerWS = evalin('caller','whos'); % arguments in caller work space
   
   % update fieldNames according to 'avoidOverWrite' flag field.
   if isfield(struct,'avoidOverWrite')
      indFieldNames = ~ismember(fieldNames,{callerWS(:).name,'avoidOverWrite'});
      fieldNames = fieldNames(indFieldNames);
      nFields = length(fieldNames);
   end
   
   if toUnpackRegular % Unpack with regular fields order
      if nargout == 0 % assign in caller
         for iField = 1:nFields
            assignin('caller',fieldNames{iField},struct.(fieldNames{iField}));
         end
      else % dump into variables
         for iField = 1:nargout
            varargout{iField} = struct.(fieldNames{iField});
         end
      end
      
   elseif toUnpackFieldNames % Unpack with fields according to fieldNames
      if nargout == 0 % assign in caller, by comparing fields to fieldNames
         for iField = 1:nFields
            assignin('caller',fieldNames{iField},struct.(fieldNames{iField}));
         end
      else % dump into variables
         assert( isequal(nFields, nargout) , ['v2struct: number of output arguments',...
            ' does not match number of field names in cell array']);
         for iField = 1:nFields
            varargout{iField} = struct.(fieldNames{iField});
         end
      end
   end
   
else  % Pack
   % issue warning for possible wrong usage when packing with an input of cell array of
   % strings without the string 'fieldNames'.
   if gotCellArrayOfStrings && ~gotFieldNames
      warning('V2STRUCT:cellArrayOfStringNotFieldNames',['v2struct.m: Please note -'...
              ' initiated packing while treating the cell array input as regular argument'])
   end
   % build cell array of input names
   callerWS = evalin('caller','whos');
   inputNames = cell(1,nargin);
   for iArgin = 1:nargin
      inputNames{iArgin} = inputname(iArgin);
   end
   nInputs = length(inputNames);
   
   % look for 'nameOfStruct2Update' variable and get the structure name
   if ~any(strcmpi(inputNames,'nameOfStruct2Update')) % no nameOfStruct2Update
      nameStructArgFound = false;
      validVarargin = varargin;
   else % nameOfStruct2Update found
      nameStructArgFound = true;
      nameStructArgLoc = strcmp(inputNames,'nameOfStruct2Update');
      nameOfStruct2Update = varargin{nameStructArgLoc};
      % valid varargin with just the inputs to pack and fieldNames if exists
      validVarargin = varargin(~strcmpi(inputNames,'nameOfStruct2Update'));
      % valid inputNames with just the inputs name to pack and fieldNames if exists
      inputNames = inputNames(~strcmpi(inputNames,'nameOfStruct2Update'));
      nInputs = length(inputNames);
      % copy structure from caller workspace to enable its updating
      if ismember(nameOfStruct2Update,{callerWS(:).name}) % verify existance
         S = evalin('caller',nameOfStruct2Update);
      else
         error(['v2struct: Bad input. Structure named ''',nameOfStruct2Update,...
            ''' was not found in workspace'])
      end
   end
   
   % fieldNames cell array exists in input
   if gotFieldNames
      nVarToPack = length(varargin)-1-double(nameStructArgFound);
      if nVarToPack == 0 % no variables to pack
         for iField = 1:nFields
            S.(fieldNames{iField}) = evalin('caller',fieldNames{iField});
         end
         
         % else - variables to pack exist
         % check for correct number of fields vs. variables to pack
      elseif ~isequal(nFields,nVarToPack)
         error(['v2struct: Bad input. Number of strings in fieldNames does not match',...
            'number of input arguments for packing.'])
      else
         for iField = 1:nFields
            S.(fieldNames{iField}) = validVarargin{iField};%varargin{iField};
         end
      end
      
   else % when input is only variables and perhaps also nameOfStruct2Update
      for iInput = 1:nInputs
         assert( ~isempty(inputNames{iInput}), ['v2struct: bad input in argument no.',...
            int2str(iArgin),'. Explicit arguments can only be called with a matching'...
            ' field names cell array']);
         S.(inputNames{iInput}) = validVarargin{iInput};
      end
   end % gotFieldNames
   
   if nargout == 0
      assignin( 'caller', 'Sv2struct',S );
   else
      varargout{1} = S;
   end
end % if nargin