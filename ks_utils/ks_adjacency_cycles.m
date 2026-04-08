% -------------------------------------------------------------------------
% SYMBOLIC PERIODIC CYCLES FROM AN 8-REGION TRANSITION GRAPH
% -------------------------------------------------------------------------
%
% DESCRIPTION
% This script defines the 8-region transition graph associated with a
% symbolic partition, constructs its adjacency matrix, and computes all
% symbolic cycles of a prescribed period.
%
% A transition A(i,j) = 1 means that the image of region R_i contains
% region R_j, that is,
%
%     F(R_i) ⊇ R_j.
%
% The script then:
%   1. Builds the adjacency matrix
%   2. Finds all symbolic cycles of a chosen period
%   3. Removes cyclic duplicates
%   4. Displays the resulting unique symbolic words
%
% AUTHOR
% Jason J. Bramburger
% -------------------------------------------------------------------------

%% Clean workspace

clear all
close all
clc

%% User parameters

nRegions = 8;       % Number of symbolic regions
targetPeriod = 4;   % Desired symbolic period

%% Build adjacency matrix

% Adjacency convention:
%   A(i,j) = 1  if  F(R_i) contains R_j
A = zeros(nRegions);

% Region-to-region transitions
A(1,[3 4])   = 1;   % R1 -> R3, R4
A(2,[1 2 3]) = 1;   % R2 -> R1, R2, R3
A(3,2)       = 1;   % R3 -> R2
A(4,[3 4 5]) = 1;   % R4 -> R3, R4, R5
A(5,[6 7])   = 1;   % R5 -> R6, R7
A(6,8)       = 1;   % R6 -> R8
A(7,[5 6 7]) = 1;   % R7 -> R5, R6, R7
A(8,[2 3 4]) = 1;   % R8 -> R2, R3, R4

%% Display adjacency matrix

disp('Adjacency matrix A:')
disp(A)

%% Find symbolic cycles of the requested period

allCycles = find_symbolic_cycles(A, targetPeriod);
uniqueCycles = remove_cyclic_duplicates(allCycles);

%% Display results

fprintf('\nAll symbolic period-%d cycles:\n', targetPeriod)
disp(allCycles)

fprintf('Unique symbolic period-%d cycles (up to cyclic rotation):\n', targetPeriod)
disp(uniqueCycles)

fprintf('Number of unique symbolic period-%d cycles: %d\n', ...
    targetPeriod, size(uniqueCycles,1))

%% Local functions

function cycles = find_symbolic_cycles(A, period)
    %FIND_SYMBOLIC_CYCLES Find all symbolic cycles of a given period.
    %
    %   cycles = find_symbolic_cycles(A, period)
    %
    % Inputs:
    %   A       - adjacency matrix
    %   period  - desired cycle length
    %
    % Output:
    %   cycles  - matrix whose rows are symbolic words of length period
    %             that close to form a cycle
    
    nStates = size(A,1);
    cycles = [];
    
        function extend_word(word)
            % Recursively extend a symbolic word until it reaches
            % the target length.
            if length(word) == period
                if A(word(end), word(1)) == 1
                    cycles = [cycles; word]; %#ok<AGROW>
                end
                return
            end
    
            currentState = word(end);
            nextStates = find(A(currentState,:));
    
            for nextState = nextStates
                extend_word([word, nextState]); %#ok<AGROW>
            end
        end
    
    for startState = 1:nStates
        extend_word(startState);
    end

end

function uniqueWords = remove_cyclic_duplicates(words)
    %REMOVE_CYCLIC_DUPLICATES Remove cyclically equivalent symbolic words.
    %
    %   uniqueWords = remove_cyclic_duplicates(words)
    %
    % Rows that differ only by cyclic rotation are treated as representing
    % the same symbolic cycle. For example,
    %
    %   [1 2 3], [2 3 1], [3 1 2]
    %
    % are identified as one cycle.
    
    if isempty(words)
        uniqueWords = words;
        return
    end
    
    [nWords, wordLength] = size(words);
    canonicalWords = zeros(nWords, wordLength);
    
    for i = 1:nWords
        word = words(i,:);
    
        % Generate all cyclic rotations
        rotations = zeros(wordLength, wordLength);
        for j = 1:wordLength
            rotations(j,:) = circshift(word, [0, -(j-1)]);
        end
    
        % Choose the lexicographically smallest rotation as canonical
        bestRotation = rotations(1,:);
        for j = 2:wordLength
            if is_lexicographically_smaller(rotations(j,:), bestRotation)
                bestRotation = rotations(j,:);
            end
        end
    
        canonicalWords(i,:) = bestRotation;
    end
    
    % Remove duplicate canonical representatives
    [~, keepIdx] = unique(canonicalWords, 'rows', 'stable');
    uniqueWords = canonicalWords(keepIdx,:);
    
    end

function tf = is_lexicographically_smaller(a, b)
    %IS_LEXICOGRAPHICALLY_SMALLER Return true if row vector a < b in lex order.
    
    tf = false;
    
    for k = 1:length(a)
        if a(k) < b(k)
            tf = true;
            return
        elseif a(k) > b(k)
            tf = false;
            return
        end
    end

end