/*
  PETPVC package for the Nix package manager.

  `default.nix` is used for building with nix-build.
  See `derivation.nix` for full description.

  Licensing
  ---------

  This file is distributed as part of PETPVC.

  Author: Ashley Gillman

  Copyright 2018 Commonwealth Scientific and Industrial Research
                 Organisation's Australian eHealth Research Centre

  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
*/

{ pkgs ? import <nixpkgs> {} }:

pkgs.callPackage ./derivation.nix { itk = pkgs.itk4; }
