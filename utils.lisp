(with-current-buffer "computeMPSATrans.m"
  (setq entries nil)
  (save-excursion
    (beginning-of-buffer)
    (while (re-search-forward (rx (or "crossTable" "SparseTensor()")) nil t)
      (push (buffer-substring-no-properties (line-beginning-position) (line-end-position)) entries)
      )
    ;; (beginning-of-buffer)
    ;; (while (re-search-forward "SparseTensor()" nil t)
      ;; (push (buffer-substring-no-properties (line-beginning-position) (line-end-position)) entries)
      ;; )
    )
  (setq entries (nreverse entries))
  )


(--map (insert it "\n") entries)

