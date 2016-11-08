function [ ] = ValidateInputs( inputs )
%ValidateInputs Validate Required inputs
%   Required inputs: 
%       * inputs.data
%       * inputs.labels
%       * inputs.rapidPTLibraryPath
%       * inputs.testingType

    assert(isfield(inputs,'data'), 'Input Error: inputs.data is a required input..');
    assert(isfield(inputs,'testingType'), 'Input Error: inputs.testingType is a required input (select between OneSample and TwoSample).');
    if(strcmp(inputs.testingType,'TwoSample'))
     %   assert(isfield(inputs,'labels'), 'Input Error: inputs.labels is a required input..');
    end

end

