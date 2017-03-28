require 'torch'
require 'nn'
require 'image'
require 'optim'

require 'loadcaffe'
require 'libcuda_utils'

require 'cutorch'
require 'cunn'

local cmd = torch.CmdLine()

-- Basic options
cmd:option('-style_image', 'examples/inputs/seated-nude.jpg', 'Style target image')
cmd:option('-content_image', 'examples/inputs/tubingen.jpg','Content target image')
cmd:option('-style_seg', '', 'Style segmentation')
cmd:option('-style_seg_idxs', '', 'Style seg idxs')
cmd:option('-content_seg', '', 'Content segmentation')
cmd:option('-content_seg_idxs', '', 'Content seg idxs')

cmd:option('-gpu', 0, 'Zero-indexed ID of the GPU to use; for CPU mode set -gpu = -1')

-- Optimization options
cmd:option('-content_weight', 5e0)
cmd:option('-style_weight', 1e2)
cmd:option('-tv_weight', 1e-3)
cmd:option('-num_iterations', 1000)

-- Output options
cmd:option('-print_iter', 1)
cmd:option('-save_iter', 100)
cmd:option('-output_image', 'out.png') 
cmd:option('-index', 1)
cmd:option('-serial', 'serial_example') 

-- Other options
cmd:option('-proto_file', 'models/VGG_ILSVRC_19_layers_deploy.prototxt')
cmd:option('-model_file', 'models/VGG_ILSVRC_19_layers.caffemodel')
cmd:option('-backend', 'nn', 'nn|cudnn|clnn')
cmd:option('-cudnn_autotune', false)
cmd:option('-seed', 612)

cmd:option('-content_layers', 'relu4_2', 'layers for content')
cmd:option('-style_layers',   'relu1_1,relu2_1,relu3_1,relu4_1,relu5_1', 'layers for style')

local function main(params)
  cutorch.setDevice(params.gpu + 1)
  cutorch.setHeapTracking(true)

  torch.manualSeed(params.seed)

  idx = cutorch.getDevice()
  print('gpu, idx = ', params.gpu, idx)

  -- content: pitie transferred input image
  local content_image = image.load(params.content_image, 3)
  local content_image_caffe = preprocess(content_image):float():cuda()
  local content_layers = params.content_layers:split(",")
 
  -- style: target model image
  local style_image = image.load(params.style_image, 3)
  local style_image_caffe = preprocess(style_image):float():cuda()
  local style_layers = params.style_layers:split(",")

  local c, h, w = content_image:size(1), content_image:size(2), content_image:size(3)
  local _, h2, w2 = style_image:size(1), style_image:size(2), style_image:size(3)
  local index = params.index

  -- segmentation images
  --[
  local content_seg = image.load(params.content_seg, 3)
  content_seg = image.scale(content_seg, w, h, 'bilinear')
  local style_seg = image.load(params.style_seg, 3)
  style_seg = image.scale(style_seg, w2, h2, 'bilinear')
  local color_codes = {'blue', 'green', 'black', 'white', 'red', 'yellow', 'grey', 'lightblue', 'purple'}
  local color_content_masks, color_style_masks = {}, {}
  for j = 1, #color_codes do
    local content_mask_j = ExtractMask(content_seg, color_codes[j])
    local style_mask_j = ExtractMask(style_seg, color_codes[j])
    table.insert(color_content_masks, content_mask_j)
    table.insert(color_style_masks, style_mask_j)
  end 
  --]]

  -- Set up the network, inserting style and content loss modules
  local content_losses, style_losses = {}, {}
  local next_content_idx, next_style_idx = 1, 1
  local net = nn.Sequential()
  
  if params.tv_weight > 0 then
    local tv_mod = nn.TVLoss(params.tv_weight):float():cuda()
    net:add(tv_mod)
  end
  
  -- load VGG-19 network
  local cnn = loadcaffe.load(params.proto_file, params.model_file, params.backend):float():cuda()

  paths.mkdir(tostring(params.serial))
  print('Exp serial:', params.serial)

  for i = 1, #cnn do
    if next_content_idx <= #content_layers or next_style_idx <= #style_layers then
      local layer = cnn:get(i)
      local name = layer.name
      local layer_type = torch.type(layer)
      local is_pooling = (layer_type == 'nn.SpatialMaxPooling' or layer_type == 'cudnn.SpatialMaxPooling')
      local is_conv    = (layer_type == 'nn.SpatialConvolution' or layer_type == 'cudnn.SpatialConvolution')
     
      net:add(layer)

      if is_pooling then
        for k = 1, #color_codes do
          color_content_masks[k] = image.scale(color_content_masks[k], math.ceil(color_content_masks[k]:size(2)/2), math.ceil(color_content_masks[k]:size(1)/2))
          color_style_masks[k]   = image.scale(color_style_masks[k],   math.ceil(color_style_masks[k]:size(2)/2),   math.ceil(color_style_masks[k]:size(1)/2))
        end
      elseif is_conv then
        local sap = nn.SpatialAveragePooling(3,3,1,1,1,1):float()
        for k = 1, #color_codes do
          color_content_masks[k] = sap:forward(color_content_masks[k]:repeatTensor(1,1,1))[1]:clone()
          color_style_masks[k]   = sap:forward(color_style_masks[k]:repeatTensor(1,1,1))[1]:clone()
        end
      end 
      color_content_masks = deepcopy(color_content_masks)
      color_style_masks = deepcopy(color_style_masks)


      if name == content_layers[next_content_idx] then
        print("Setting up content layer", i, ":", layer.name)
        local target = net:forward(content_image_caffe):clone()
        local loss_module = nn.ContentLoss(params.content_weight, target, false):float():cuda()
        net:add(loss_module)
        table.insert(content_losses, loss_module)
        next_content_idx = next_content_idx + 1
      end

     if name == style_layers[next_style_idx] then
        print("Setting up style layer  ", i, ":", layer.name)
        local gram = GramMatrix():float():cuda()
        local target_features = net:forward(style_image_caffe):clone()

        local target_grams = {}

        for j = 1, #color_codes do 
          local l_style_mask_ori = color_style_masks[j]:clone():cuda()
          local l_style_mask = l_style_mask_ori:repeatTensor(1,1,1):expandAs(target_features)
          local l_style_mean = l_style_mask_ori:mean()
           
          local masked_target_features = torch.cmul(l_style_mask, target_features)
          local masked_target_gram = gram:forward(masked_target_features):clone()
          if l_style_mean > 0 then
            masked_target_gram:div(target_features:nElement() * l_style_mean)
          end 
          table.insert(target_grams, masked_target_gram)
        end 
        
        local loss_module = nn.StyleLossWithSeg(params.style_weight, target_grams, color_content_masks, color_codes, next_style_idx, false):float():cuda()

        net:add(loss_module)
        table.insert(style_losses, loss_module)
        next_style_idx = next_style_idx + 1
      end 

    end
  end

  -- We don't need the base CNN anymore, so clean it up to save memory.
  cnn = nil
  for i=1,#net.modules do
    local module = net.modules[i]
    if torch.type(module) == 'nn.SpatialConvolutionMM' then
        -- remove these, not used, but uses gpu memory
        module.gradWeight = nil
        module.gradBias = nil
    end
  end
  collectgarbage()

  local mean_pixel = torch.CudaTensor({103.939, 116.779, 123.68})
  local meanImage = mean_pixel:view(3, 1, 1):expandAs(content_image_caffe)

  local img = torch.randn(content_image:size()):float():mul(0.0001):cuda()

  -- Run it through the network once to get the proper size for the gradient
  -- All the gradients will come from the extra loss modules, so we just pass
  -- zeros into the top of the net on the backward pass.
  local y = net:forward(img)
  local dy = img.new(#y):zero()

  -- Declaring this here lets us access it in maybe_print
  local optim_state = {
      maxIter = params.num_iterations,
      tolX = 0, tolFun = -1,
      verbose=true, 
  }

  local function maybe_print(t, loss)
    local verbose = (params.print_iter > 0 and t % params.print_iter == 0)
    if verbose then
      print(string.format('Iteration %d / %d', t, params.num_iterations))
      for i, loss_module in ipairs(content_losses) do
        print(string.format('  Content %d loss: %f', i, loss_module.loss))
      end
      for i, loss_module in ipairs(style_losses) do
        print(string.format('  Style %d loss: %f', i, loss_module.loss))
      end
      print(string.format('  Total loss: %f', loss))
    end
  end

  local function maybe_save(t)
    local should_save = params.save_iter > 0 and t % params.save_iter == 0
    should_save = should_save or t == params.num_iterations
    if should_save then
      local disp = deprocess(img:double())
      disp = image.minmax{tensor=disp, min=0, max=1}
      local filename = params.serial .. '/out' .. tostring(index) .. '_t_' .. tostring(t) .. '.png'
      image.save(filename, disp)
    end
  end

  local num_calls = 0
  local function feval(AffineModel) 
    num_calls = num_calls + 1

    local output = torch.add(img, meanImage)
    local input  = torch.add(content_image_caffe, meanImage)

    net:forward(img)
   
    local gradient_VggNetwork = net:updateGradInput(img, dy)

    local grad = gradient_VggNetwork

    local loss = 0
    for _, mod in ipairs(content_losses) do
      loss = loss + mod.loss
    end
    for _, mod in ipairs(style_losses) do
      loss = loss + mod.loss
    end
    maybe_print(num_calls, loss)
    maybe_save(num_calls)
    
    collectgarbage()

    -- optim.lbfgs expects a vector for gradients
    return loss, grad:view(grad:nElement()) 
  end
  
  -- Run optimization.
  local x, losses = optim.lbfgs(feval, img, optim_state)  
end
 

function build_filename(output_image, iteration)
  local ext = paths.extname(output_image)
  local basename = paths.basename(output_image, ext)
  local directory = paths.dirname(output_image)
  return string.format('%s/%s_%d.%s',directory, basename, iteration, ext)
end

-- Preprocess an image before passing it to a Caffe model.
-- We need to rescale from [0, 1] to [0, 255], convert from RGB to BGR,
-- and subtract the mean pixel.
function preprocess(img)
  local mean_pixel = torch.DoubleTensor({103.939, 116.779, 123.68})
  local perm = torch.LongTensor{3, 2, 1}
  img = img:index(1, perm):mul(256.0)
  mean_pixel = mean_pixel:view(3, 1, 1):expandAs(img)
  img:add(-1, mean_pixel)
  return img
end


-- Undo the above preprocessing.
function deprocess(img)
  local mean_pixel = torch.DoubleTensor({103.939, 116.779, 123.68})
  mean_pixel = mean_pixel:view(3, 1, 1):expandAs(img)
  img = img + mean_pixel
  local perm = torch.LongTensor{3, 2, 1}
  img = img:index(1, perm):div(256.0)
  return img
end

function deepcopy(orig)
    local orig_type = type(orig)
    local copy
    if orig_type == 'table' then
        copy = {}
        for orig_key, orig_value in next, orig, nil do
            copy[deepcopy(orig_key)] = deepcopy(orig_value)
        end
        setmetatable(copy, deepcopy(getmetatable(orig)))
    else -- number, string, boolean, etc
        copy = orig
    end
    return copy
end

-- Define an nn Module to compute content loss in-place
local ContentLoss, parent = torch.class('nn.ContentLoss', 'nn.Module')

function ContentLoss:__init(strength, target, normalize)
  parent.__init(self)
  self.strength = strength
  self.target = target
  self.normalize = normalize or false
  self.loss = 0
  self.crit = nn.MSECriterion()
end

function ContentLoss:updateOutput(input)
  if input:nElement() == self.target:nElement() then
    self.loss = self.crit:forward(input, self.target) * self.strength
  else
    print('WARNING: Skipping content loss')
  end
  self.output = input
  return self.output
end

function ContentLoss:updateGradInput(input, gradOutput)
  if input:nElement() == self.target:nElement() then
    self.gradInput = self.crit:backward(input, self.target)
  end
  if self.normalize then
    self.gradInput:div(torch.norm(self.gradInput, 1) + 1e-8)
  end
  self.gradInput:mul(self.strength)
  self.gradInput:add(gradOutput)
  return self.gradInput
end

-- Returns a network that computes the CxC Gram matrix from inputs
-- of size C x H x W
function GramMatrix()
  local net = nn.Sequential()
  net:add(nn.View(-1):setNumInputDims(2))
  local concat = nn.ConcatTable()
  concat:add(nn.Identity())
  concat:add(nn.Identity())
  net:add(concat)
  net:add(nn.MM(false, true))
  return net
end


-- Define an nn Module to compute style loss in-place
local StyleLoss, parent = torch.class('nn.StyleLoss', 'nn.Module')

function StyleLoss:__init(strength, target, normalize)
  parent.__init(self)
  self.normalize = normalize or false
  self.strength = strength
  self.target = target
  self.loss = 0
  
  self.gram = GramMatrix()
  self.G = nil
  self.crit = nn.MSECriterion()
end

function StyleLoss:updateOutput(input)
  self.G = self.gram:forward(input)
  self.G:div(input:nElement())
  self.loss = self.crit:forward(self.G, self.target)
  self.loss = self.loss * self.strength
  self.output = input
  return self.output
end

function StyleLoss:updateGradInput(input, gradOutput)
  local dG = self.crit:backward(self.G, self.target)
  dG:div(input:nElement())
  self.gradInput = self.gram:backward(input, dG)
  if self.normalize then
    self.gradInput:div(torch.norm(self.gradInput, 1) + 1e-8)
  end
  self.gradInput:mul(self.strength)
  self.gradInput:add(gradOutput)
  return self.gradInput
end


function ExtractMask(seg, color)
  local mask = nil
  if color == 'green' then 
    mask = torch.lt(seg[1], 0.1)
    mask:cmul(torch.gt(seg[2], 1-0.1))
    mask:cmul(torch.lt(seg[3], 0.1))
  elseif color == 'black' then 
    mask = torch.lt(seg[1], 0.1)
    mask:cmul(torch.lt(seg[2], 0.1))
    mask:cmul(torch.lt(seg[3], 0.1))
  elseif color == 'white' then
    mask = torch.gt(seg[1], 1-0.1)
    mask:cmul(torch.gt(seg[2], 1-0.1))
    mask:cmul(torch.gt(seg[3], 1-0.1))
  elseif color == 'red' then 
    mask = torch.gt(seg[1], 1-0.1)
    mask:cmul(torch.lt(seg[2], 0.1))
    mask:cmul(torch.lt(seg[3], 0.1))
  elseif color == 'blue' then
    mask = torch.lt(seg[1], 0.1)
    mask:cmul(torch.lt(seg[2], 0.1))
    mask:cmul(torch.gt(seg[3], 1-0.1))
  elseif color == 'yellow' then
    mask = torch.gt(seg[1], 1-0.1)
    mask:cmul(torch.gt(seg[2], 1-0.1))
    mask:cmul(torch.lt(seg[3], 0.1))
  elseif color == 'grey' then 
    mask = torch.cmul(torch.gt(seg[1], 0.5-0.1), torch.lt(seg[1], 0.5+0.1))
    mask:cmul(torch.cmul(torch.gt(seg[2], 0.5-0.1), torch.lt(seg[2], 0.5+0.1)))
    mask:cmul(torch.cmul(torch.gt(seg[3], 0.5-0.1), torch.lt(seg[3], 0.5+0.1)))
  elseif color == 'lightblue' then
    mask = torch.lt(seg[1], 0.1)
    mask:cmul(torch.gt(seg[2], 1-0.1))
    mask:cmul(torch.gt(seg[3], 1-0.1))
  elseif color == 'purple' then 
    mask = torch.gt(seg[1], 1-0.1)
    mask:cmul(torch.lt(seg[2], 0.1))
    mask:cmul(torch.gt(seg[3], 1-0.1))
  else 
    print('ExtractMask(): color not recognized, color = ', color)
  end 
  return mask:float()
end

-- Define style loss with segmentation 
local StyleLossWithSeg, parent = torch.class('nn.StyleLossWithSeg', 'nn.Module')

--function StyleLossWithSeg:__init(strength, target_grams, color_content_masks, content_seg_idxs, layer_id, normalize)
function StyleLossWithSeg:__init(strength, target_grams, color_content_masks, color_codes, layer_id, normalize)
  parent.__init(self)
  self.strength = strength
  self.target_grams = target_grams
  self.color_content_masks = deepcopy(color_content_masks)
  self.color_codes = color_codes
  --self.content_seg_idxs = content_seg_idxs
  self.normalize = normalize

  self.loss = 0
  self.gram = GramMatrix()
  self.crit = nn.MSECriterion()

  self.layer_id = layer_id
end 

function StyleLossWithSeg:updateOutput(input)
  self.output = input
  return self.output
end 

function StyleLossWithSeg:updateGradInput(input, gradOutput)
  self.loss = 0
  self.gradInput = gradOutput:clone()
  self.gradInput:zero()
  for j = 1, #self.color_codes do 
    local l_content_mask_ori = self.color_content_masks[j]:clone():cuda()
    local l_content_mask = l_content_mask_ori:repeatTensor(1,1,1):expandAs(input) 
    local l_content_mean = l_content_mask_ori:mean()
    
    local masked_input_features = torch.cmul(l_content_mask, input)
    local masked_input_gram = self.gram:forward(masked_input_features):clone()
    if l_content_mean > 0 then 
      masked_input_gram:div(input:nElement() * l_content_mean)
    end

    local loss_j = self.crit:forward(masked_input_gram, self.target_grams[j])
    loss_j = loss_j * self.strength * l_content_mean
    self.loss = self.loss + loss_j

    local dG = self.crit:backward(masked_input_gram, self.target_grams[j])
 
    dG:div(input:nElement())

    local gradient = self.gram:backward(masked_input_features, dG) 

    if self.normalize then 
      gradient:div(torch.norm(gradient, 1) + 1e-8)
    end

    self.gradInput:add(gradient)
  end   

  self.gradInput:mul(self.strength)
  self.gradInput:add(gradOutput)
  return self.gradInput
end 


local TVLoss, parent = torch.class('nn.TVLoss', 'nn.Module')

function TVLoss:__init(strength)
  parent.__init(self)
  self.strength = strength
  self.x_diff = torch.Tensor()
  self.y_diff = torch.Tensor()
end

function TVLoss:updateOutput(input)
  self.output = input
  return self.output
end

-- TV loss backward pass inspired by kaishengtai/neuralart
function TVLoss:updateGradInput(input, gradOutput)
  self.gradInput:resizeAs(input):zero()
  local C, H, W = input:size(1), input:size(2), input:size(3)
  self.x_diff:resize(3, H - 1, W - 1)
  self.y_diff:resize(3, H - 1, W - 1)
  self.x_diff:copy(input[{{}, {1, -2}, {1, -2}}])
  self.x_diff:add(-1, input[{{}, {1, -2}, {2, -1}}])
  self.y_diff:copy(input[{{}, {1, -2}, {1, -2}}])
  self.y_diff:add(-1, input[{{}, {2, -1}, {1, -2}}])
  self.gradInput[{{}, {1, -2}, {1, -2}}]:add(self.x_diff):add(self.y_diff)
  self.gradInput[{{}, {1, -2}, {2, -1}}]:add(-1, self.x_diff)
  self.gradInput[{{}, {2, -1}, {1, -2}}]:add(-1, self.y_diff)
  self.gradInput:mul(self.strength)
  self.gradInput:add(gradOutput)
  return self.gradInput
end

function TVGradient(input, gradOutput, strength)
  local C, H, W = input:size(1), input:size(2), input:size(3)
  local gradInput = torch.CudaTensor(C, H, W):zero()
  local x_diff = torch.CudaTensor()
  local y_diff = torch.CudaTensor()
  x_diff:resize(3, H - 1, W - 1)
  y_diff:resize(3, H - 1, W - 1)
  x_diff:copy(input[{{}, {1, -2}, {1, -2}}])
  x_diff:add(-1, input[{{}, {1, -2}, {2, -1}}])
  y_diff:copy(input[{{}, {1, -2}, {1, -2}}])
  y_diff:add(-1, input[{{}, {2, -1}, {1, -2}}])
  gradInput[{{}, {1, -2}, {1, -2}}]:add(x_diff):add(y_diff)
  gradInput[{{}, {1, -2}, {2, -1}}]:add(-1, x_diff)
  gradInput[{{}, {2, -1}, {1, -2}}]:add(-1, y_diff)
  gradInput:mul(strength)
  gradInput:add(gradOutput)
  return gradInput
end 
 
local params = cmd:parse(arg)
main(params)
    
