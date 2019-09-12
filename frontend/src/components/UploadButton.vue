<template>
  <span>
    <v-btn @click="trigger" color="primary" elevation="1" small>
      <v-icon small left>mdi-cloud-upload</v-icon>
      Upload
    </v-btn>
    <input :multiple="multiple" class="visually-hidden" type="file" v-on:change="files" ref="fileInput" />
  </span>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { Component, Emit, Prop, Vue } from "vue-property-decorator";

@Component
export default class UploadButton extends Vue {
  experimentContext = experimentModule.context(this.$store);

  @Prop(Number) id!: number;
  @Prop({ default: false }) multiple!: boolean;

  @Emit()
  async files(e): Promise<FileList> {
    const formData = new FormData();
    const file = e.target.files[0];
    formData.append("file", file, file.name);
    e.target.value = "";
    await this.experimentContext.actions.upload({ id: this.id, data: formData });
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
