<template>
  <div class="editor">
    <div class="container" style="height: 90vh">
      <div id="editor" class="node-editor" ></div>
    </div>
    <div class="dock"></div>
  </div>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { settingsModule } from '@/modules/settings';
  import { NumComponent } from '@/views/main/experiment/workflow/NumComponent';
  import Rete from 'rete';
  import CommentPlugin from 'rete-comment-plugin';
  import ConnectionPlugin from 'rete-connection-plugin';
  import ContextMenuPlugin from 'rete-context-menu-plugin';
  import DockPlugin from 'rete-dock-plugin';
  import MinimapPlugin from 'rete-minimap-plugin';
  import VueRenderPlugin from 'rete-vue-render-plugin';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class WorkflowView extends Vue {
    readonly experimentContext = experimentModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    mounted() {
      this.showDiagram();
    }

    private showDiagram() {
      const container = document.getElementById('editor');
      if (!container) {
        return;
      }

      const editor = new Rete.NodeEditor('demo@0.1.0', container);

      editor.use(DockPlugin, {
        container: document.getElementById('editor'),
        itemClass: 'item',
        plugins: [VueRenderPlugin]
      });

      editor.use(ConnectionPlugin);
      editor.use(VueRenderPlugin);
      // editor.use(DockPlugin);
      editor.use(MinimapPlugin);
      editor.use(ContextMenuPlugin, {
        searchBar: false, // true by default
        searchKeep: title => true, // leave item when searching, optional. For example, title => ['Refresh'].includes(title)
        delay: 100,
        allocate(component) {
          return ['Submenu'];
        },
        rename(component) {
          return component.name;
        },
        items: {
          'Click me'() {
            console.log('Works!');
          },
        },
        nodeItems: {
          'Click me'() {
            console.log('Works for node!');
          },
        },
      });
      editor.use(CommentPlugin, {
        margin: 20, // default indent for new frames is 30px
      });

      const numComponent = new NumComponent();
      editor.register(numComponent);

      const engine = new Rete.Engine('demo@0.1.0');
      engine.register(numComponent);

      editor.on(['process', 'nodecreated', 'noderemoved', 'connectioncreated', 'connectionremoved'], async () => {
        await engine.abort();
        await engine.process(editor.toJSON());
      });
    }
  }
</script>

<style scoped>
  .editor {
    display: flex;
    flex-wrap: nowrap;
    flex-direction: column;
    height: 100vh;
  }

  .dock {
    height: 100px;
    overflow-x: auto;
    overflow-y: hidden;
    white-space: nowrap;
  }

  .dock-item {
    display: inline-block;
    vertical-align: top;
    transform: scale(0.8);
    transform-origin: 50% 0;
  }

  .container {
    flex: 1;
    overflow: hidden;
  }
</style>
