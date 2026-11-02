function signature = capsule_source_signature(root)
    % CAPSULE_SOURCE_SIGNATURE Hash MATLAB sources before reusing any checkpoint.
    %   A change to a source file, including its comments, starts a different
    %   run. This keeps checkpoints tied to their producing source files.

    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    files = dir(fullfile(root, '**', '*.m'));
    paths = sort(string(fullfile({files.folder}, {files.name})));
    digest = java.security.MessageDigest.getInstance('SHA-256');
    for k = 1:numel(paths)
        relative = erase(paths(k), string(root) + filesep);
        digest.update(typecast(unicode2native(char(relative), 'UTF-8'), 'int8'));
        file = fopen(paths(k), 'rb');
        if file < 0
            error('MDDIPIC:ReadSource', 'Cannot read %s.', paths(k));
        end
        guard = onCleanup(@() fclose(file));
        digest.update(typecast(fread(file, Inf, '*uint8'), 'int8'));
        clear guard;
    end
    signature = string(lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, [])));
end
