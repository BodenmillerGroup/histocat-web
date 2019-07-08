<template>
  <v-tooltip top>
    <span>Upload Slide</span>
    <v-btn slot="activator" flat @click="trigger">
      <v-icon>mdi-cloud-upload</v-icon>
    </v-btn>
    <input :multiple="multiple" class="visually-hidden" type="file" v-on:change="files" ref="fileInput">
  </v-tooltip>
</template>

<script lang="ts">
  import { experimentModule } from '@/modules/experiment';
  import { Component, Emit, Prop, Vue } from 'vue-property-decorator';

  @Component
  export default class UploadButton extends Vue {
    experimentContext = experimentModule.context(this.$store);

    @Prop(Number) id!: number;
    @Prop(String) color!: string;
    @Prop({ default: false }) multiple!: boolean;

    @Emit()
    async files(e): Promise<FileList> {
      const formData = new FormData();
      const file = e.target.files[0];
      formData.append('file', file, file.name);
      await this.experimentContext.actions.uploadSlide({ id: this.id, data: formData });
      return e.target.files;
    }

    trigger() {
      (this.$refs.fileInput as HTMLElement).click();
    }
  }
</script>

<style scoped>
  .visually-hidden {
    position: absolute !important;
    height: 1px;
    width: 1px;
    overflow: hidden;
    clip: rect(1px, 1px, 1px, 1px);
  }
</style>
