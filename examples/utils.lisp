(defun matlab-to-python-syntax (beg end)
  (interactive "r")
  (save-excursion
    (goto-char beg)
    (save-excursion (while (search-forward "**" end t)
                      (replace-match ".^")))
    (save-excursion (while (search-forward "*" end t)
                      (replace-match ".*")))
    )
  )

