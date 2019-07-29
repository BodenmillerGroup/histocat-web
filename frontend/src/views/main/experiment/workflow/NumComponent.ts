import Rete, { Node } from 'rete';
import { NodeData, WorkerInputs, WorkerOutputs } from 'rete/types/core/data';

const numSocket = new Rete.Socket('Number value');

export class NumComponent extends Rete.Component {
  constructor() {
    super('Number');
  }

  builder(node: Node): any {
    const out = new Rete.Output('num', 'Number', numSocket);
    node.addOutput(out);
  }

  worker(node: NodeData, inputs: WorkerInputs, outputs: WorkerOutputs) {
    outputs['num'] = node.data.num;
  }
}
