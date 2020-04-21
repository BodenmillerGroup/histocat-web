import rambda from "rambda";

const DATA_KEYS = [
  "class",
  "staticClass",
  "style",
  "attrs",
  "props",
  "domProps",
  "on",
  "nativeOn",
  "directives",
  "scopesSlots",
  "slot",
  "ref",
  "key",
];

function mutateKey(key) {
  return "" + key + `-cloned-cid`;
}

function extractData(vnode, isComp) {
  const data = rambda.pick(DATA_KEYS, vnode.data);
  if (isComp) {
    const cOpts = vnode.componentOptions;
    Object.assign(data, {
      props: cOpts.propsData,
      on: cOpts.listeners,
    });
  }

  if (data.key) {
    data.key = mutateKey(data.key);
  }

  return data;
}

export function cloneElement(vnode, newData = {}) {
  // use the context that the original vnode was created in.
  const h = vnode.context && vnode.context.$createElement;
  const isComp = !!vnode.componentOptions;
  const isText = !vnode.tag; // this will also match comments but those will be dropped, essentially
  const children = isComp ? vnode.componentOptions.children : vnode.children;

  if (isText) return vnode.text;

  const data = extractData(vnode, isComp);

  const tag = isComp ? vnode.componentOptions.Ctor : vnode.tag;

  const childNodes = children ? children.map((c) => cloneElement(c)) : undefined;
  return h(tag, data, childNodes);
}
