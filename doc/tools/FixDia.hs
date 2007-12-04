module Main where

import Text.Regex
import System

subs = lines >>> split "\n" >>> 

main = do
  args  <- getArgs
  lines <- readFile $ args !! 1
  pairs <- (split "\n" >>> )
  interact $ mkRegex